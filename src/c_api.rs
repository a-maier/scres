use std::os::raw::{c_double, c_void};

use cres::{
    c_api::{
        cres::Search,
        event::{EventView, TypeSetView},
    },
    distance::EuclWithScaledPt,
    event::{Event, EventBuilder},
    traits::Distance,
};
use noisy_float::prelude::*;
use particle_id::ParticleID;

use crate::resampler::{Resampler, ResamplerBuilder};

/// Resampling options
#[repr(C)]
#[derive(Copy, Clone, Debug)]
pub struct Opt {
    /// Nearest-neighbour search algorithm
    neighbour_search: Search,
    /// Extra contribution to distance proportional to difference in pt
    ///
    /// This parameter is ignored when using a custom distance. Otherwise,
    /// it corresponds to the τ parameter of
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
        .distance(dist.clone())
        .neighbour_search(neighbour_search)
        .build();
    let resampler: &mut dyn CResampler = Box::leak(Box::new(resampler));
    Box::into_raw(Box::new(resampler)) as _
}

/// Delete a resampler
#[no_mangle]
pub unsafe extern "C" fn scres_free(scres: *mut c_void) {
    assert!(!scres.is_null());
    let _ = Box::from_raw(scres as *mut &mut dyn CResampler);
}

/// Reserve space for events (optional)
#[no_mangle]
pub unsafe extern "C" fn scres_reserve(scres: *mut c_void, cap: usize) {
    let scres = scres as *mut &mut dyn CResampler;
    (*scres).reserve(cap);
}

/// Add an event
#[no_mangle]
pub unsafe extern "C" fn scres_push_event(scres: *mut c_void, event: EventView) {
    let scres = scres as *mut &mut dyn CResampler;
    (*scres).push(event);
}

/// Remove the *last* event and retrieve its weights
///
/// Returns a null pointer if there are no events left.
#[no_mangle]
#[must_use]
pub unsafe extern "C" fn scres_next_weights(scres: *mut c_void) -> *const c_double {
    let scres = scres as *mut &mut dyn CResampler;
    match (*scres).next_weights() {
        Some(wts) => wts.as_ptr(),
        None => std::ptr::null(),
    }
}

/// Construct a cell with the `n`th event as seed and resample
#[no_mangle]
pub unsafe extern "C" fn scres_resample(
    scres: *const c_void,
    seed: usize,
    max_cell_size: c_double
) {
    let scres = scres as *const &mut dyn CResampler;
    (*scres).resample_cell(seed, max_cell_size);
}

pub trait CResampler {
    fn resample_cell(&self, seed: usize, max_cell_size: f64);

    fn reserve(&mut self, cap: usize);

    fn push(&mut self, event: EventView);

    fn next_weights(&mut self) -> Option<&[f64]>;
}

impl<D: Distance + Send + Sync> CResampler for Resampler<D> {
    fn resample_cell(&self, seed: usize, max_cell_size: f64) {
        self.resample_cell(seed, n64(max_cell_size));
    }

    fn reserve(&mut self, cap: usize) {
        self.reserve(cap);
    }

    fn push(&mut self, event: EventView) {
        self.push(ToEvent(event).into());
    }

    fn next_weights(&mut self) -> Option<&[f64]> {
        self.next_weights()
            .map(|w| unsafe { std::mem::transmute(w) })
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
        // TODO: reserve
        let mut event = EventBuilder::new();
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
