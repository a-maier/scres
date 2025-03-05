use cres::{
    c_api::cres::Search,
    cell::Cell,
    distance::{DistWrapper, EuclWithScaledPt},
    event::Event,
    neighbour_search::{NaiveNeighbourSearch, TreeSearch},
    traits::{Distance, NeighbourSearchAlgo},
    N64,
};

/// A single-cell resampler
#[derive(Debug)]
pub struct Resampler<D> {
    distance: D,
    neighbour_search: Search,
    events: Vec<Event>,
    last_retrieved_weights: Vec<N64>,
}

/// Resample the cell with the event number `seed` as cell seed
impl<D> Resampler<D>
where
    D: Distance + Send + Sync,
{
    pub fn resample_cell(&self, seed: usize, max_cell_size: N64) {
        let mut cell = match self.neighbour_search {
            Search::Tree => {
                let neighbour_search = TreeSearch::new_with_dist(
                    self.events.len(),
                    DistWrapper::new(&self.distance, &self.events),
                    max_cell_size,
                );
                Cell::new(&self.events, seed, &neighbour_search)
            }
            Search::Naive => {
                let neighbour_search = NaiveNeighbourSearch::new_with_dist(
                    self.events.len(),
                    DistWrapper::new(&self.distance, &self.events),
                    max_cell_size,
                );
                Cell::new(&self.events, seed, &neighbour_search)
            }
        };
        cell.resample();
    }
}

impl<D> Resampler<D> {
    /// Reserve space for `cap` events
    pub fn reserve(&mut self, cap: usize) {
        self.events.reserve(cap);
    }

    /// Add an event
    pub fn push(&mut self, event: Event) {
        self.events.push(event)
    }

    /// Retrieve the weights of the given event
    pub fn get_weights(&mut self, pos: usize) -> &[N64] {
        self.last_retrieved_weights =
            self.events[pos].weights.read().iter().copied().collect();
        &self.last_retrieved_weights
    }

    /// Retrieve the number of weights of the given event
    pub fn get_num_weights(&self, pos: usize) -> usize {
        self.events[pos].n_weights()
    }

    /// Retrieve the weights of the given event
    pub fn set_weights(&self, pos: usize, weights: &[N64]) {
        // TODO: this is awkward
        // the natural way would be `copy_from_slice`,
        // but the `Weights` struct does not expose it
        let mut target_weights = self.events[pos].weights.write();
        for (lhs, rhs) in target_weights.iter_mut().zip(weights) {
            *lhs = *rhs
        }
    }

    /// Remove all events
    pub fn clear(&mut self) {
        self.events.clear();
    }
}

/// Construct a [Resampler] object
pub struct ResamplerBuilder<D> {
    distance: D,
    neighbour_search: Search,
}

impl Default for ResamplerBuilder<EuclWithScaledPt> {
    fn default() -> Self {
        Self {
            distance: Default::default(),
            neighbour_search: Search::Tree,
        }
    }
}

impl<D> ResamplerBuilder<D> {
    /// Set the nearest neighbour search algorithm
    pub fn neighbour_search(
        self,
        neighbour_search: Search,
    ) -> ResamplerBuilder<D> {
        ResamplerBuilder {
            distance: self.distance,
            neighbour_search,
        }
    }

    /// Set the distance
    pub fn distance<DD>(self, distance: DD) -> ResamplerBuilder<DD> {
        ResamplerBuilder {
            distance,
            neighbour_search: self.neighbour_search,
        }
    }

    /// Build a [Resampler]
    pub fn build(self) -> Resampler<D> {
        let Self {
            distance,
            neighbour_search,
        } = self;
        Resampler {
            distance,
            neighbour_search,
            events: vec![],
            last_retrieved_weights: vec![],
        }
    }
}
