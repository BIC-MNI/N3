#include "FitField.h"

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "../SplineSmooth/fieldIO.h"

/* Declared in splineSmooth.cc, which has no header of its own for them. */
DblMat volume_domain(VIO_Volume volume);
DblMat reduced_domain(VIO_Volume volume);
void fitSplinesToVolume(Spline *spline, VIO_Volume volume, int subsample);
void fitSplinesToVolume(Spline *spline, VIO_Volume volume, VIO_Volume mask_volume,
                        const DblMat &domain, int subsample);
void fitSplinesToVolumeLookup(TBSplineVolume *spline, VIO_Volume volume,
                              int subsample);
void fitSplinesToVolumeLookup(TBSplineVolume *spline, VIO_Volume volume,
                              VIO_Volume mask_volume, const DblMat &domain,
                              int subsample);

namespace n3 {

Field::Field(Spline *s, const DblMat &d, double dist, enum spline_type t)
  : spline(s), domain(d), distance(dist), type(t) {}

Field::~Field() { delete spline; }

Field *fit_field(VIO_Volume volume, VIO_Volume mask, enum spline_type type,
                 double distance, double lambda, int subsample)
{
  /* The driver passes -full_support with -b_spline and omits it with
   * -tp_spline (nu_estimate_np_and_em.in:571-579), which is what sends the two
   * bases down splineSmooth.cc:124-130's two branches. */
  bool full_support = (type == b_spline);
  DblMat domain = (!full_support && mask) ? reduced_domain(mask)
                                          : volume_domain(volume);

  Spline *spline;
  if(type == b_spline)
    {
      VIO_Real separations[VIO_N_DIMENSIONS];
      int sizes[VIO_N_DIMENSIONS];
      VIO_Real start[VIO_N_DIMENSIONS] = { 0.0, 0.0, 0.0 };
      get_volume_separations(volume, separations);
      get_volume_sizes(volume, sizes);

      TBSplineVolume *b = new TBSplineVolume(domain, start, separations, sizes,
                                             distance, lambda);
      if(mask) fitSplinesToVolumeLookup(b, volume, mask, domain, subsample);
      else     fitSplinesToVolumeLookup(b, volume, subsample);
      spline = b;
    }
  else
    {
      spline = createThinPlateSpline(domain, distance, lambda, FALSE);
      if(mask) fitSplinesToVolume(spline, volume, mask, domain, subsample);
      else     fitSplinesToVolume(spline, volume, subsample);
    }

  return new Field(spline, domain, distance, type);
}

void evaluate_field(Field *field, VIO_Volume target, VIO_Volume mask)
{
  double real_min, real_max;

  if(field->type == b_spline)
    {
      /* The basis has to be laid out on the target's grid: TBSplineVolume
       * caches per-grid lookup tables.  The coefficients carry across because
       * the number of them depends on the domain and the knot spacing alone,
       * and the spline's coordinates are index times step from voxel (0,0,0).
       * This is what inputCompactField does (fieldIO.cc:312-326), with the
       * coefficients coming from memory instead of from a file. */
      VIO_Real separations[VIO_N_DIMENSIONS];
      int sizes[VIO_N_DIMENSIONS];
      VIO_Real start[VIO_N_DIMENSIONS] = { 0.0, 0.0, 0.0 };
      get_volume_separations(target, separations);
      get_volume_sizes(target, sizes);

      TBSplineVolume *on_target =
        new TBSplineVolume(field->domain, start, separations, sizes,
                           field->distance, 1.0, FALSE);
      if(on_target->putCoefficients(field->spline->getCoefficients()) == FALSE)
        {
          fprintf(stderr, "n3::evaluate_field: coefficient count mismatch\n");
          exit(1);
        }

      if(mask)
        smoothVolumeLookup(on_target, target, mask, &real_min, &real_max);
      else
        smoothVolumeLookup(on_target, target, &real_min, &real_max);

      delete on_target;
    }
  else
    {
      /* A thin plate spline is evaluated at a point, so it needs no per-grid
       * basis and carries across unchanged. */
      if(mask)
        smoothVolume(field->spline, target, mask, &real_min, &real_max);
      else
        smoothVolume(field->spline, target, &real_min, &real_max);
    }
}

void save_field(const std::string &path, Field *field, VIO_Volume like,
                const std::string &command)
{
  if(outputCompactField(path.c_str(), field->domain, field->distance,
                        field->spline->getCoefficients(), field->type,
                        command.c_str(), like) != VIO_OK)
    {
      fprintf(stderr, "n3::save_field: cannot write %s\n", path.c_str());
      exit(1);
    }
}

void evaluate_saved_field(const std::string &path, VIO_Volume target,
                          VIO_Volume mask)
{
  Spline *spline = NULL;
  enum spline_type type;

  if(inputCompactField((char *) path.c_str(), &spline, &type, target) != VIO_OK)
    {
      fprintf(stderr, "n3::evaluate_saved_field: cannot read %s\n", path.c_str());
      exit(1);
    }

  double real_min, real_max;
  if(type == b_spline)
    {
      if(mask)
        smoothVolumeLookup((TBSplineVolume *) spline, target, mask,
                           &real_min, &real_max);
      else
        smoothVolumeLookup((TBSplineVolume *) spline, target,
                           &real_min, &real_max);
    }
  else
    {
      if(mask) smoothVolume(spline, target, mask, &real_min, &real_max);
      else     smoothVolume(spline, target, &real_min, &real_max);
    }

  delete spline;
}

/* --- .imp reading --------------------------------------------------------
 * The MNI field file format is four keyword lines then a row of coefficients;
 * fieldIO.cc::inputCompactField is the reader, but it keeps Distance and
 * Domain to itself.  The three helpers below pull just those two fields out
 * of the same text so a Field can be reconstructed on another grid.  They
 * deliberately do not parse the coefficients -- inputCompactField does that. */

/* The number that ends the Distance = n; line. */
static double imp_distance(const char *path)
{
  FILE *f = fopen(path, "r");
  if(!f) { fprintf(stderr, "n3::load_field: cannot open %s\n", path); exit(1); }
  double distance = 0.0;
  char line[512];
  while(fgets(line, sizeof(line), f))
    if(sscanf(line, "Distance = %lf", &distance) == 1) break;
  fclose(f);
  return distance;
}

/* The six world coordinates that make up the Domain block, stopping at the
 * semicolon.  They may spread over several lines, which is why they are read
 * with strtod across the buffer rather than with sscanf per line. */
static void imp_domain_world(const char *path, VIO_Real world[6])
{
  FILE *f = fopen(path, "r");
  if(!f) { fprintf(stderr, "n3::load_field: cannot open %s\n", path); exit(1); }

  char line[1024];
  while(fgets(line, sizeof(line), f))
    if(strstr(line, "Domain"))
      {
        break;
      }

  int got = 0;
  while(got < 6 && fgets(line, sizeof(line), f))
    {
      char *p = line;
      char *end;
      while(got < 6 && *p && *p != ';')
        {
          while(*p && (isspace((unsigned char) *p) || *p == ';')) p++;
          if(!*p || *p == ';') break;
          world[got++] = strtod(p, &end);
          if(end == p) break;
          p = end;
        }
    }
  fclose(f);
  if(got != 6)
    {
      fprintf(stderr, "n3::load_field: malformed Domain in %s\n", path);
      exit(1);
    }
}

Field *load_field(const std::string &path, VIO_Volume like)
{
  Spline *spline = NULL;
  enum spline_type type;

  if(inputCompactField((char *) path.c_str(), &spline, &type, like) != VIO_OK)
    {
      fprintf(stderr, "n3::load_field: cannot read %s\n", path.c_str());
      exit(1);
    }

  double distance = imp_distance(path.c_str());

  VIO_Real world[6];
  imp_domain_world(path.c_str(), world);

  /* inputCompactField converts the world coordinates back to voxels and
   * scales by the separations, giving a domain in mm on the volume's grid
   * (:285-292).  Replicated here because Field needs that domain to lay the
   * basis out on the target grid (evaluate_field). */
  /* Field file Domain blocks write one axis per line, x0 x1, y0 y1, z0 z1
   * (outputCompactField, fieldIO.cc:124-132), so the six reals are
   * interleaved: [x0, x1, y0, y1, z0, z1].  inputCompactField reads them the
   * same way (fieldIO.cc:259-268). */
  VIO_Real world0[VIO_N_DIMENSIONS] = { world[0], world[2], world[4] };
  VIO_Real world1[VIO_N_DIMENSIONS] = { world[1], world[3], world[5] };

  VIO_Real separations[VIO_N_DIMENSIONS];
  VIO_Real voxel0[VIO_N_DIMENSIONS], voxel1[VIO_N_DIMENSIONS];
  get_volume_separations(like, separations);
  convert_world_to_voxel(like, world0[0], world0[1], world0[2], voxel0);
  convert_world_to_voxel(like, world1[0], world1[1], world1[2], voxel1);

  DblMat domain(VIO_N_DIMENSIONS, 2);
  for(int i = 0; i < VIO_N_DIMENSIONS; i++)
    {
      domain(i, 0) = voxel0[i] * separations[i];
      domain(i, 1) = voxel1[i] * separations[i];
    }

  return new Field(spline, domain, distance, type);
}

}  // namespace n3
