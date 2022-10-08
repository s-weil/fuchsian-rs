use std::ops::Deref;

use crate::{
    algebraic_extensions::{Group, MulIdentity},
    set_extensions::{SetRestriction, Wrapper},
};

pub trait Determinant<T> {
    fn det(&self) -> T;
}

impl<M, S, T> Determinant<T> for M
where
    M: Wrapper<Inner = S>,
    S: Determinant<T>,
{
    fn det(&self) -> T {
        self.deref().det()
    }
}

/// We identify linear transformations satisfying the condition `determinant == 1` with the
/// subset [SL(2,R)](https://en.wikipedia.org/wiki/SL2(R)) within 2x2 matrices.
/// Due to mathematical limitations in rust we cannot model this condition,
/// hence we check the condition upon construction (SetRestriction) and may assume the condition in the following.
pub trait SpecialLinear<T>: Determinant<T> + SetRestriction {
    /// overwrite condition to satisfy `determinant == 1`
    fn condition(&self) -> bool
    where
        T: MulIdentity + PartialEq,
    {
        self.det() == T::one()
    }
    // TODO: check if overwrite actually works
}

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
                        item = item.inv();
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

impl<Space, G, IG> Action<Space> for G
where
    G: Deref<Target = IG>,
    IG: Action<Space>,
{
    fn map(&self, x: &Space) -> Space {
        self.deref().map(x)
    }
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
