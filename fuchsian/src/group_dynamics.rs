use crate::{
    algebraic_extensions::{Group, MulIdentity},
    set_extensions::SetRestriction,
};

pub trait Determinant<T> {
    fn determinant(&self) -> T;
}

pub trait SpecialLinear<T>: Determinant<T> + SetRestriction {
    fn restriction(s: &Self) -> bool
    where
        T: PartialEq,
        Self: MulIdentity,
    {
        let one: Self = MulIdentity::one();
        s.determinant() == one.determinant()
    }
}

// use crate::algebraic_extensions::AddIdentity;
// // The [mathematical group](https://en.wikipedia.org/wiki/Group_(mathematics)#Definition)
// // definition except for the associativity identity.
// // /// In particular, `MoebiusTransformation<T>` is an additive group
// /// ```
// /// use std::ops::{Add, Neg};
// /// impl<T> Group for T
// /// where
// ///     T: AddIdentityElement + Add<Output = Self> + Neg<Output = Self> + Copy + PartialEq + Sized,
// /// {
// ///     fn combine(&self, other: &Self) -> Self {
// ///         *self + *other
// ///     }
// ///
// ///     fn identity() -> Self {
// ///         Self::zero()
// ///     }
// ///
// ///     fn inverse(&self) -> Self {
// ///         -*self
// ///     }
// /// }
// /// ```
// pub trait Group: PartialEq + Sized {
//     /// The binary operation.
//     fn combine(&self, other: &Self) -> Self;

//     /// The identity element.
//     fn identity() -> Self;

//     /// The inverse element.
//     fn inverse(&self) -> Self;

//     // how to model identities in general?
//     fn associativity_check(&self, b: &Self, c: &Self) -> bool {
//         let ab_c = self.combine(b).combine(c);
//         let a_bc = self.combine(&b.combine(c));
//         ab_c == a_bc
//     }
// }

/// A [finitely generated group](https://en.wikipedia.org/wiki/Finitely_generated_group) `G` has
/// some finite generating set so that every `g \in G`can be written as a finite combination
/// of these elements (and their inveserses).
pub trait FinitelyGeneratedGroup {
    type GroupElement: Group;

    fn generators(&self) -> &[Self::GroupElement];
    // fn inverse_generators(&self) -> &[Self::GroupElement];

    // TODO: add inverse set?
    // fn inverse_generators(&self) -> &[Self::GroupElement] {
    //     self.generators().iter()
    // }
}

#[derive(Default)]
pub enum PickGeneratorMode {
    #[default]
    Sequential,
    Random,
    // Single(a) ?
}

struct GeneratorPicker<'a, G>
where
    G: FinitelyGeneratedGroup,
{
    mode: PickGeneratorMode,
    group: &'a G,
    cursor: usize,
    max_n: usize,
    generator_len: usize,
}

impl<'a, G> GeneratorPicker<'a, G>
where
    G: FinitelyGeneratedGroup,
    G::GroupElement: Clone,
{
    fn new(mode: PickGeneratorMode, group: &'a G, max_n: usize) -> Self {
        Self {
            mode,
            max_n,
            cursor: 0,
            generator_len: group.generators().len(),
            group,
        }
    }
}

impl<'a, G> Iterator for GeneratorPicker<'a, G>
where
    G: FinitelyGeneratedGroup,
    G::GroupElement: Clone,
{
    type Item = G::GroupElement;

    fn next(&mut self) -> Option<Self::Item> {
        match self.mode {
            PickGeneratorMode::Sequential => {
                if self.cursor < self.max_n {
                    let mut item: Self::Item =
                        self.group.generators()[(self.cursor % self.generator_len)].clone();

                    // TODO: precalculate inverses!
                    if self.cursor % (2 * self.generator_len) >= self.generator_len {
                        item = item.inverse();
                    }

                    self.cursor += 1;
                    return Some(item);
                }
                None
            }
            PickGeneratorMode::Random => {
                todo!("create random nr within 0..length and choose random element")
            }
        }
    }
}

pub trait Action<Space> {
    fn map(&self, x: &Space) -> Space;
}

/// The mathematical (left) [group action](https://en.wikipedia.org/wiki/Group_action) of
/// a group `G` acting on a space `X`.
pub trait GroupAction<Space>: Action<Space> + Group
where
    Space: PartialEq,
{
    // TODO: how to model identities in general?
    fn identity_check(&self, x: &Space) -> bool {
        (&Group::identity() as &Self).map(x).eq(x)
    }

    fn compatibility_check(&self, h: &Self, x: &Space) -> bool {
        let g_hx = self.map(&h.map(x));
        let gh_x = (self.combine(h)).map(x);
        g_hx.eq(&gh_x)
    }
}

// TODO: could instead implement iterator for Orbit
pub struct Orbit<Space>
where
    Space: Sized,
{
    pub points: Vec<Space>,
    // generators: &'a G,
}

impl<Space> Orbit<Space>
where
    Space: Sized,
{
    pub fn create<Group>(
        group: &Group,
        base_point: &Space,
        n_points: usize,
        pick_generator: Option<PickGeneratorMode>,
    ) -> Self
    where
        Group: FinitelyGeneratedGroup,
        Group::GroupElement: Action<Space> + Clone,
        Space: Clone,
    {
        let mut points = Vec::with_capacity(n_points);
        let pick_mode = pick_generator.unwrap_or_default();

        let generator = GeneratorPicker::new(pick_mode, group, n_points);
        let mut point_cursor = base_point.clone();

        for g in generator {
            point_cursor = g.map(&point_cursor);
            points.push(point_cursor.clone());
        }

        Self {
            points,
            // generators: group,
        }
    }
}

// TODO:
// - revisit trait definitions and required bounds
// - Fundamental Domain
// - check performance of group action / orbit. i.e. what is faster, combining group elements and then action, or iterative actions
