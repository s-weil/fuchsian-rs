use rand::{
    distributions::{DistIter, Uniform},
    rngs::ThreadRng,
};

use crate::{
    algebraic_extensions::{Group, MulIdentity},
    set_extensions::{SetRestriction, Wrapper},
};
use std::ops::Deref;

pub trait Determinant<T> {
    fn det(&self) -> T;
}

impl<W, I, T> Determinant<T> for W
where
    W: Wrapper<Inner = I>,
    I: Determinant<T>,
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

#[derive(Default, PartialEq, Eq)]
pub enum PickGeneratorMode {
    #[default]
    Sequential,
    Random,
    // Single(a) ?
}

struct GeneratorPicker<'a, G> {
    mode: PickGeneratorMode,
    generators: &'a Vec<G>,
    cursor: usize,
    max_n: usize,
    generator_len: usize,
    // rand_iter: Option<DistIter<Uniform<usize>, ThreadRng, usize>>,
    rand_iter: DistIter<Uniform<usize>, ThreadRng, usize>,
}

impl<'a, G> GeneratorPicker<'a, G>
where
    G: Clone,
{
    fn new(mode: PickGeneratorMode, generators: &'a Vec<G>, max_n: usize) -> Self {
        let generator_len = generators.len();
        let rand_iter = random_iter(generator_len);

        Self {
            mode,
            max_n,
            cursor: 0,
            generator_len,
            generators,
            rand_iter,
        }
    }
}

// TODO: make it feature dependent
fn random_iter(u_bound: usize) -> DistIter<Uniform<usize>, ThreadRng, usize> {
    use rand::{thread_rng, Rng};
    let rng = thread_rng();
    let gen_range = Uniform::new(0, u_bound);
    rng.sample_iter(gen_range)
}

impl<'a, G> Iterator for GeneratorPicker<'a, G>
where
    G: Clone,
{
    type Item = G;

    fn next(&mut self) -> Option<Self::Item> {
        match self.mode {
            PickGeneratorMode::Sequential => {
                if self.cursor < self.max_n {
                    let item = self.generators[(self.cursor % self.generator_len)].clone();
                    self.cursor += 1;
                    return Some(item);
                }
                None
            }
            PickGeneratorMode::Random => {
                if self.cursor < self.max_n {
                    // let grp_idx = self.rand_iter.as_mut().unwrap().next().unwrap();
                    let grp_idx = self.rand_iter.next().unwrap();
                    let item = self.generators[grp_idx].clone();
                    self.cursor += 1;
                    return Some(item);
                }
                None
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
    pub fn sample<Group>(
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

        let mut generators = Vec::with_capacity(2 * group.generators().len());
        for g in group.generators().iter() {
            generators.push(g.clone());
            generators.push(g.inv().clone());
        }

        let generator = GeneratorPicker::new(pick_mode, &generators, n_points);
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
