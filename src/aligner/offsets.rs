use std::fmt::Debug;
use std::hash::Hash;
use std::ops::{BitAnd, Not, Shr};

use num::{FromPrimitive, Unsigned, One, Bounded};
use num::traits::{SaturatingAdd, SaturatingSub};

pub trait OffsetType: FromPrimitive + Unsigned + PartialEq + Eq
    + PartialOrd + Ord + Default + Copy + Hash + Debug + Bounded + SaturatingSub + SaturatingAdd
    + Not<Output=Self> + BitAnd<Output=Self> + Shr<Output=Self>
{
    fn new(value: usize) -> Self;
    fn as_usize(&self) -> usize;
    fn as_isize(&self) -> isize;
    fn increase_one(&self) -> Self;
}

impl OffsetType for u8 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}

impl OffsetType for u16 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }

    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}

impl OffsetType for u32 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }
    
    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}

impl OffsetType for u64 {
    #[inline(always)]
    fn new(value: usize) -> Self {
        value as Self
    }

    #[inline(always)]
    fn as_usize(&self) -> usize {
        *self as usize
    }
    
    #[inline(always)]
    fn as_isize(&self) -> isize {
        *self as isize
    }

    #[inline(always)]
    fn increase_one(&self) -> Self {
        *self + Self::one()
    }
}

#[cfg(test)]
mod tests {
    use super::OffsetType;

    fn check_offset<T: OffsetType + PartialEq + From<u8>>(value: u8) {
        let off: T = <T as OffsetType>::new(value as usize);
        assert_eq!(off, T::from(value));
        assert_eq!(off.as_usize(), value as usize);
        assert_eq!(off.as_isize(), value as isize);
        assert_eq!(off.increase_one(), T::from(value + 1));
    }

    #[test]
    fn u8_offsets() {
        check_offset::<u8>(8);
    }

    #[test]
    fn u16_offsets() {
        check_offset::<u16>(16);
    }

    #[test]
    fn u32_offsets() {
        check_offset::<u32>(32);
    }
}
