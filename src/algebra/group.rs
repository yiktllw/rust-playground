pub trait Semigroup {
    type Item;
    fn combine(a: &Self::Item, b: &Self::Item) -> Self::Item;
}

pub trait Monoid: Semigroup {
    fn identity() -> Self::Item;
}

pub trait Group: Monoid {
    fn inverse(a: &Self::Item) -> Self::Item;
}

pub trait AbelianGroup: Group {
    fn commutative(a: &Self::Item, b: &Self::Item) -> bool;
}