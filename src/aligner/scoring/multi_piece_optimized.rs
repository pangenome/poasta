use super::{Score, AlignState};
use rustc_hash::FxHashMap;

/// Optimized visited cell for multi-piece gap states
/// Uses inline storage for common cases to avoid HashMap overhead
#[derive(Clone, Debug)]
pub struct OptimizedMultiPieceCell {
    // Basic states
    visited_m: Score,
    visited_i: Score,
    visited_d: Score,
    visited_i2: Score,
    visited_d2: Score,
    
    // Inline storage for multi-piece states (up to 4 piece/position combinations)
    // This covers most practical cases without HashMap overhead
    inline_states: [(u16, Score); 4],
    inline_count: u8,
    
    // Overflow storage for rare cases with many states
    overflow: Option<Box<FxHashMap<u16, Score>>>,
}

impl Default for OptimizedMultiPieceCell {
    fn default() -> Self {
        Self {
            visited_m: Score::Unvisited,
            visited_i: Score::Unvisited,
            visited_d: Score::Unvisited,
            visited_i2: Score::Unvisited,
            visited_d2: Score::Unvisited,
            inline_states: [(0, Score::Unvisited); 4],
            inline_count: 0,
            overflow: None,
        }
    }
}

impl OptimizedMultiPieceCell {
    #[inline]
    pub fn get_score(&self, state: AlignState) -> Score {
        match state {
            AlignState::Match => self.visited_m,
            AlignState::Insertion => self.visited_i,
            AlignState::Deletion => self.visited_d,
            AlignState::Insertion2 => self.visited_i2,
            AlignState::Deletion2 => self.visited_d2,
            AlignState::MultiInsertion { piece, position } | 
            AlignState::MultiDeletion { piece, position } => {
                let key = ((piece as u16) << 8) | (position as u16);
                
                // First check inline storage
                for i in 0..self.inline_count as usize {
                    if self.inline_states[i].0 == key {
                        return self.inline_states[i].1;
                    }
                }
                
                // Then check overflow if it exists
                if let Some(ref overflow) = self.overflow {
                    overflow.get(&key).copied().unwrap_or(Score::Unvisited)
                } else {
                    Score::Unvisited
                }
            }
        }
    }
    
    #[inline]
    pub fn set_score(&mut self, state: AlignState, score: Score) {
        match state {
            AlignState::Match => self.visited_m = score,
            AlignState::Insertion => self.visited_i = score,
            AlignState::Deletion => self.visited_d = score,
            AlignState::Insertion2 => self.visited_i2 = score,
            AlignState::Deletion2 => self.visited_d2 = score,
            AlignState::MultiInsertion { piece, position } | 
            AlignState::MultiDeletion { piece, position } => {
                let key = ((piece as u16) << 8) | (position as u16);
                self.set_multi_score(key, score);
            }
        }
    }
    
    #[inline]
    fn set_multi_score(&mut self, key: u16, score: Score) {
        // First try to update existing entry in inline storage
        for i in 0..self.inline_count as usize {
            if self.inline_states[i].0 == key {
                self.inline_states[i].1 = score;
                return;
            }
        }
        
        // Try to add to inline storage if there's space
        if (self.inline_count as usize) < self.inline_states.len() {
            self.inline_states[self.inline_count as usize] = (key, score);
            self.inline_count += 1;
            return;
        }
        
        // Otherwise use overflow storage
        if self.overflow.is_none() {
            self.overflow = Some(Box::new(FxHashMap::default()));
        }
        self.overflow.as_mut().unwrap().insert(key, score);
    }
    
    #[inline]
    pub fn update_score_if_lower(&mut self, state: AlignState, new_score: Score) -> bool {
        let current = self.get_score(state);
        if new_score < current {
            self.set_score(state, new_score);
            true
        } else {
            false
        }
    }
}

/// Alternative: Use a more compact representation that doesn't track individual positions
/// This trades some optimality for much better performance
#[derive(Clone, Debug, Default)]
pub struct CompactMultiPieceCell {
    // Basic states
    visited_m: Score,
    visited_i: Score,
    visited_d: Score,
    
    // Track best score for each piece (not position within piece)
    // Most multi-piece models have 2-3 pieces, so this is very efficient
    visited_i_pieces: [Score; 4],
    visited_d_pieces: [Score; 4],
}

impl CompactMultiPieceCell {
    #[inline]
    pub fn get_score_compact(&self, state: AlignState) -> Score {
        match state {
            AlignState::Match => self.visited_m,
            AlignState::Insertion => self.visited_i,
            AlignState::Deletion => self.visited_d,
            AlignState::Insertion2 => self.visited_i_pieces[1],
            AlignState::Deletion2 => self.visited_d_pieces[1],
            AlignState::MultiInsertion { piece, .. } => {
                if piece < 4 {
                    self.visited_i_pieces[piece as usize]
                } else {
                    Score::Unvisited
                }
            },
            AlignState::MultiDeletion { piece, .. } => {
                if piece < 4 {
                    self.visited_d_pieces[piece as usize]
                } else {
                    Score::Unvisited
                }
            }
        }
    }
    
    #[inline]
    pub fn update_score_if_lower_compact(&mut self, state: AlignState, new_score: Score) -> bool {
        match state {
            AlignState::Match => {
                if new_score < self.visited_m {
                    self.visited_m = new_score;
                    true
                } else {
                    false
                }
            },
            AlignState::Insertion => {
                if new_score < self.visited_i {
                    self.visited_i = new_score;
                    true
                } else {
                    false
                }
            },
            AlignState::Deletion => {
                if new_score < self.visited_d {
                    self.visited_d = new_score;
                    true
                } else {
                    false
                }
            },
            AlignState::Insertion2 => {
                if new_score < self.visited_i_pieces[1] {
                    self.visited_i_pieces[1] = new_score;
                    true
                } else {
                    false
                }
            },
            AlignState::Deletion2 => {
                if new_score < self.visited_d_pieces[1] {
                    self.visited_d_pieces[1] = new_score;
                    true
                } else {
                    false
                }
            },
            AlignState::MultiInsertion { piece, .. } => {
                if piece < 4 && new_score < self.visited_i_pieces[piece as usize] {
                    self.visited_i_pieces[piece as usize] = new_score;
                    true
                } else {
                    false
                }
            },
            AlignState::MultiDeletion { piece, .. } => {
                if piece < 4 && new_score < self.visited_d_pieces[piece as usize] {
                    self.visited_d_pieces[piece as usize] = new_score;
                    true
                } else {
                    false
                }
            }
        }
    }
}