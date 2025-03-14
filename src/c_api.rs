use std::{
    os::raw::{c_double, c_void},
    ptr::slice_from_raw_parts,
};

use cres::{
    c_api::{
        cres::Search,
        event::{EventView, TypeSetView},
    },
    distance::EuclWithScaledPt,
    event::{Event, EventBuilder},
    N64, n64,
    traits::Distance,
    ParticleID,
};

use crate::resampler::{Resampler, ResamplerBuilder};

/// Resampling options
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct Opt {
    /// Nearest-neighbour search algorithm
    neighbour_search: Search,
    /// Extra contribution to distance proportional to difference in pt
    ///
    /// This parameter corresponds to the τ parameter of
    /// [arXiv:2109.07851](https://arxiv.org/abs/2109.07851)
    pt_weight: c_double,
}

/// Create a new resampler
#[no_mangle]
#[must_use]
pub extern "C" fn scres_new(opt: Opt) -> *mut c_void {
    let Opt {
        neighbour_search,
        pt_weight,
    } = opt;
    let dist = EuclWithScaledPt::new(n64(pt_weight));
    let resampler = ResamplerBuilder::default()
        .distance(dist)
        .neighbour_search(neighbour_search)
        .build();
    let resampler: &mut dyn CResampler = Box::leak(Box::new(resampler));
    Box::into_raw(Box::new(resampler)) as _
}

/// Delete a resampler
///
/// # Safety
/// The resampler must have been previous constructed with `scres_new`.
///
#[no_mangle]
pub unsafe extern "C" fn scres_free(scres: *mut c_void) {
    assert!(!scres.is_null());
    let _ = Box::from_raw(scres as *mut &mut dyn CResampler);
}

/// Reserve space for events (optional)
///
/// # Safety
/// The resampler must have been previous constructed with `scres_new`.
///
#[no_mangle]
pub unsafe extern "C" fn scres_reserve(scres: *mut c_void, cap: usize) {
    let scres = scres as *mut &mut dyn CResampler;
    (*scres).reserve(cap);
}

/// Add an event
///
/// # Safety
/// The resampler must have been previous constructed with `scres_new`.
///
#[no_mangle]
pub unsafe extern "C" fn scres_push_event(
    scres: *mut c_void,
    event: EventView,
) {
    let scres = scres as *mut &mut dyn CResampler;
    (*scres).push(event);
}

/// Construct a cell with the `n`th event as seed and resample
///
/// # Safety
/// The resampler must have been previous constructed with `scres_new`.
///
#[no_mangle]
pub unsafe extern "C" fn scres_resample(
    scres: *const c_void,
    seed: usize,
    max_cell_size: c_double,
) {
    let scres = scres as *const &mut dyn CResampler;
    (*scres).resample_cell(seed, max_cell_size);
}

/// Returns the weights of the chosen event
///
/// # Safety
/// The resampler must have been previous constructed with `scres_new`.
///
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn scres_get_weights(
    scres: *mut c_void,
    pos: usize,
) -> *const c_double {
    let scres = scres as *mut &mut dyn CResampler;
    (*scres).get_weights(pos).as_ptr()
}

/// Get the number of weights of the chosen event
///
/// # Safety
/// - The resampler must have been previous constructed with `scres_new`.
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn scres_get_num_weights(
    scres: *const c_void,
    pos: usize,
) -> usize {
    let scres = scres as *const &dyn CResampler;
    (*scres).get_num_weights(pos)
}

/// Sets the weights of the chosen event
///
/// # Safety
/// - The resampler must have been previous constructed with `scres_new`.
/// - The number of elements in the `weights` array has to be at least
///   as large as the existing number of weights in the event. Extra
///   elements will be ignored. `weights` must not be null.
#[no_mangle]
pub unsafe extern "C" fn scres_set_weights(
    scres: *const c_void,
    pos: usize,
    weights: *const c_double,
) {
    let num_weights = scres_get_num_weights(scres, pos);
    let scres = scres as *const &dyn CResampler;
    let weights = slice_from_raw_parts(weights, num_weights);
    (*scres).set_weights(pos, weights.as_ref().unwrap())
}

/// Delete all pushed events
///
/// # Safety
/// The resampler must have been previous constructed with `scres_new`.
///
#[no_mangle]
pub unsafe extern "C" fn scres_clear(scres: *mut c_void) {
    let scres = scres as *mut &mut dyn CResampler;
    (*scres).clear()
}

pub trait CResampler {
    fn resample_cell(&self, seed: usize, max_cell_size: f64);

    fn reserve(&mut self, cap: usize);

    fn push(&mut self, event: EventView);

    fn get_weights(&mut self, pos: usize) -> &[f64];

    fn get_num_weights(&self, pos: usize) -> usize;

    fn set_weights(&self, pos: usize, weights: &[f64]);

    fn clear(&mut self);
}

impl<D: Distance + Send + Sync> CResampler for Resampler<D> {
    fn resample_cell(&self, seed: usize, max_cell_size: f64) {
        self.resample_cell(seed, n64(max_cell_size));
    }

    fn reserve(&mut self, cap: usize) {
        self.reserve(cap);
    }

    fn push(&mut self, event: EventView) {
        // TODO: remove `ToEvent` once cres 0.9 is available
        self.push(ToEvent(event).into());
    }

    fn get_num_weights(&self, pos: usize) -> usize {
        Resampler::get_num_weights(self, pos)
    }

    fn get_weights(&mut self, pos: usize) -> &[f64] {
        // Safety: N64 and f64 have the same memory layout and alignment
        unsafe { std::mem::transmute(self.get_weights(pos)) }
    }

    fn set_weights(&self, pos: usize, weights: &[f64]) {
        // Safety: N64 and f64 have the same memory layout and alignment
        let weights: &[N64] = unsafe { std::mem::transmute(weights) };
        self.set_weights(pos, weights)
    }

    fn clear(&mut self) {
        self.clear();
    }
}

struct ToEvent<'a>(EventView<'a>);

impl<'a> From<ToEvent<'a>> for Event {
    fn from(view: ToEvent<'a>) -> Self {
        let EventView {
            id: _,
            weights,
            type_sets,
            n_weights,
            n_type_sets,
        } = view.0;
        let n_particles = (0..n_type_sets)
            .map(|n| unsafe { (*type_sets.add(n)).n_momenta })
            .sum();
        let mut event = EventBuilder::with_capacity(n_particles);
        for n_weight in 0..n_weights {
            let weight = unsafe { *weights.add(n_weight) };
            event.add_weight(n64(weight));
        }
        for n_set in 0..n_type_sets {
            let TypeSetView {
                pid,
                momenta,
                n_momenta,
                ..
            } = unsafe { *type_sets.add(n_set) };
            for n_p in 0..n_momenta {
                let p = unsafe { *momenta.add(n_p) };
                event.add_outgoing(ParticleID::new(pid), p.map(n64).into());
            }
        }
        event.build()
    }
}

#[cfg(test)]
mod tests {
    use core::f64;

    use cres::c_api::event::TypeSet;

    use super::*;

    #[test]
    fn c_api() {
        unsafe {
            let opt = Opt {
                neighbour_search: Search::Tree,
                pt_weight: 0.0,
            };
            let resampler = scres_new(opt);
            scres_reserve(resampler, 2);

            let jets = TypeSet {
                pid: 90,
                momenta: vec![
                    [
                        0.86042412975E+02,
                        0.18299527188E+02,
                        0.50776693328E+02,
                        -0.67008593105E+02,
                    ],
                    [
                        0.80026513931E+03,
                        -0.18299527188E+02,
                        -0.50776693328E+02,
                        -0.79844295220E+03,
                    ],
                ],
            };
            let view = jets.view();
            let weights = -1.0;
            let event = EventView {
                id: 0,
                weights: &weights as _,
                type_sets: &view as _,
                n_weights: 1,
                n_type_sets: 1,
            };
            scres_push_event(resampler, event);

            let jets = TypeSet {
                pid: 90,
                momenta: vec![
                    [
                        0.49452408437E+02,
                        0.20789583719E+02,
                        -0.23718791628E+02,
                        0.38088749425E+02,
                    ],
                    [
                        0.10452662667E+03,
                        -0.20789583719E+02,
                        0.23718791628E+02,
                        0.99654542370E+02,
                    ],
                ],
            };
            let view = jets.view();
            let weights = 1.0;
            let event = EventView {
                id: 0,
                weights: &weights as _,
                type_sets: &view as _,
                n_weights: 1,
                n_type_sets: 1,
            };
            scres_push_event(resampler, event);

            scres_resample(resampler, 0, f64::MAX);

            assert_eq!(*scres_get_weights(resampler, 0), 0.0);
            assert_eq!(*scres_get_weights(resampler, 1), 0.0);

            scres_clear(resampler);

            scres_free(resampler);
        }
    }
}
