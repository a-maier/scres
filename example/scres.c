/* Example for running cell resampling on two dijet events
 *
 * To run this example, first download the scres library from
 * https://github.com/a-maier/scres/releases. Then copy the compiled
 * libraries (`libscres.a` and `libscres.so` on linux) and the
 * generated headers `scres.h` and `cres.h` to a directory where they
 * can be found by your C compiler.
 *
 * Now compile the example, for example with
 * ```
 * cc -o scres examples/cres.c -lscres -lcres -lm
 * ```
 *
 * Finally, run with
 * ```
 * ./scres
 * ```
 */
#include <assert.h>
#include <float.h>
#include <stdint.h>

#include "cres.h"
#include "scres.h"

typedef struct {
  double weight;

  FourMomentum jet_momenta[2];
} DiJetEvent;

int main() {
  /* Define the events */

  DiJetEvent events[2];
  events[0] = (DiJetEvent){
    .weight = -1.0,
    .jet_momenta = {
      {0.86042412975E+02, 0.18299527188E+02,  0.50776693328E+02, -0.67008593105E+02},
      {0.80026513931E+03, -0.18299527188E+02, -0.50776693328E+02, -0.79844295220E+03},
    }
  };
  events[1] = (DiJetEvent){
    .weight = 1.0,
    .jet_momenta = {
      {0.49452408437E+02, 0.20789583719E+02, -0.23718791628E+02,  0.38088749425E+02},
      {0.10452662667E+03, -0.20789583719E+02, 0.23718791628E+02, 0.99654542370E+02}
    }
  };

  /* Create the resampler */
  ScresOpt opt = (ScresOpt) {
    .neighbour_search = Tree,
    .pt_weight = 0.0
  };
  void * resampler = scres_new(opt);

  /* reserve space for two events (optional) */
  scres_reserve(resampler, 2u);

  /* add first event */
  TypeSetView jet_view = (TypeSetView) {
    /* Particle id
     *
     * Ultimately, scres does not care about it, as long as different
     * particle types have different ids
     */
    .pid = 90,
    .n_momenta = 2u, /* number of particles of this type, aka. jets */
    .momenta = events[0].jet_momenta /* particle momenta */
  };
  EventView event = (EventView) {
    .id = 0, /* this field is ignored */
    .n_type_sets = 1u, /* number of different particle types */
    .n_weights = 1u, /* number of weights */
    .weights = &events[0].weight /* event weights */,
  };
  scres_push_event(resampler, event);

  /* add second event */
  jet_view.momenta = events[1].jet_momenta;
  event.weights = &events[1].weight;

  /* resample with first event as seed and unlimited cell size */
  scres_resample(resampler, 0, DBL_MAX);

  /* delete events in *reverse* order, retrieving their weights */
  double weight = *scres_next_weights(resampler);
  assert(weight == 0.0);
  weight = *scres_next_weights(resampler);
  assert(weight == 0.0);

  /* no more events left */
  assert(scres_next_weights(resampler) == NULL);

  /* clean up */
  scres_free(resampler);
}
