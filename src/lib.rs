#![warn(missing_docs)]

//! The `stabby` crate provides an implementation of the data structure
//! and associated algorithms presented in:
//! 
//! > Schmidt, Jens M. "Interval stabbing problems in small integer ranges."
//! > In Algorithms and Computation: 20th International Symposium, ISAAC 2009,
//! > Honolulu, Hawaii, USA, December 16-18, 2009. Proceedings 20, pp. 163-172. 
//! > Springer Berlin Heidelberg, 2009.
//! 
//! The data is built from a sorted list of [`Interval`](crate::Stabby) objects
//! each representing a closed interval over unsigned integers. Once constructed,
//! the data structure supports efficient queries.
//! 
//! The implementation uses two layers. The outer layer maps the sparse domain of
//! intervals over `u64` into a dense domain using a rank/select data structure,
//! then the intervals over the dense domain are used to build the data structre
//! described by Schmidt.
//! 
//! # Examples
//! 
//! ```rust
//! use stabby::Interval;
//! use stabby::Stabby;

//! let mutyh = vec![
//!     Interval::new(45_329_163, 45_329_437), 
//!     Interval::new(45_330_516, 45_330_557),
//!     Interval::new(45_331_182, 45_331_334),
//!     Interval::new(45_331_420, 45_331_556),
//!     Interval::new(45_332_763, 45_332_834),
//!     Interval::new(45_334_391, 45_334_511),
//!     Interval::new(45_340_219, 45_340_447),
//! ];
//! let idx = Stabby::new(&mutyh);
//! assert!(idx.stabs(45_331_258));
//! assert!(!idx.stabs(45_332_840));
//! assert_eq!(idx.stab(45_331_258),
//!            vec![Interval::new(45_331_182, 45_331_334)]);
//! assert_eq!(idx.stab_interval(&Interval::new(45_331_151, 45_331_880)),
//!            vec![Interval::new(45_331_182, 45_331_334),
//!                 Interval::new(45_331_420, 45_331_556)]);
//! ```

mod listy;
mod dense;
mod sparse;

pub use sparse::Interval;
pub use sparse::Stabby;
