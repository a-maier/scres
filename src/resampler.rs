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

    /// Remove the *last* event and retrieve its weights
    pub fn next_weights(&mut self) -> Option<&[N64]> {
        let next_ev = self.events.pop()?;
        self.last_retrieved_weights =
            next_ev.weights.into_inner().iter().copied().collect();
        Some(&self.last_retrieved_weights)
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
