use crate::{
    algebraic_extensions::{Group, MulIdentity},
    set_extensions::{SetRestriction, Wrapper},
};
use rand::{
    distributions::{DistIter, Uniform},
    rngs::ThreadRng,
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
}

#[derive(Default, PartialEq, Eq)]
pub enum PickGeneratorMode {
    #[default]
    Sequential,
    Random,
}

struct SequentialPicker<'a, G> {
    generators: &'a Vec<G>,
    cursor: usize,
    max_n: usize,
    generator_len: usize,
}

impl<'a, G> SequentialPicker<'a, G>
where
    G: Clone,
{
    fn new(generators: &'a Vec<G>, max_n: usize) -> Self {
        let generator_len = generators.len();

        Self {
            max_n,
            cursor: 0,
            generator_len,
            generators,
        }
    }
}

impl<'a, G> Iterator for SequentialPicker<'a, G>
where
    G: Clone,
{
    type Item = G;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cursor < self.max_n {
            let item = self.generators[(self.cursor % self.generator_len)].clone();
            self.cursor += 1;
            return Some(item);
        }
        None
    }
}

struct RandomPicker<'a, G> {
    generators: &'a Vec<G>,
    cursor: usize,
    max_n: usize,
    rand_iter: DistIter<Uniform<usize>, ThreadRng, usize>,
}

impl<'a, G> RandomPicker<'a, G>
where
    G: Clone,
{
    fn new(generators: &'a Vec<G>, max_n: usize) -> Self {
        let rand_iter = random_iter(generators.len());

        Self {
            max_n,
            cursor: 0,
            rand_iter,
            generators,
        }
    }
}

impl<'a, G> Iterator for RandomPicker<'a, G>
where
    G: Clone,
{
    type Item = G;

    fn next(&mut self) -> Option<Self::Item> {
        if self.cursor < self.max_n {
            let grp_idx = self.rand_iter.next().unwrap();
            let item = self.generators[grp_idx].clone();
            self.cursor += 1;
            return Some(item);
        }
        None
    }
}

// TODO: make it feature dependent
fn random_iter(u_bound: usize) -> DistIter<Uniform<usize>, ThreadRng, usize> {
    use rand::{thread_rng, Rng};
    let rng = thread_rng();
    let gen_range = Uniform::new(0, u_bound);
    rng.sample_iter(gen_range)
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

        let mut generators = Vec::with_capacity(2 * group.generators().len());
        if group.generators().len() > 1 {
            // order of adding generators is important, so that an element with its inverse don't cancel each other immediately
            for g in group.generators().iter() {
                generators.push(g.clone());
            }
            for g in group.generators().iter() {
                generators.push(g.inv().clone());
            }
        } else if let Some(g) = group.generators().iter().next() {
            generators.push(g.clone())
        }

        let mut point_cursor = base_point.clone();

        match pick_generator.unwrap_or_default() {
            PickGeneratorMode::Sequential => {
                let generator = SequentialPicker::new(&generators, n_points);
                for g in generator {
                    point_cursor = g.map(&point_cursor);
                    points.push(point_cursor.clone());
                }
            }
            PickGeneratorMode::Random => {
                let generator = RandomPicker::new(&generators, n_points);
                for g in generator {
                    point_cursor = g.map(&point_cursor);
                    points.push(point_cursor.clone());
                }
            }
        };

        Self { points }
    }
}

#[cfg(test)]
mod tests {
    use super::SequentialPicker;

    #[test]
    fn test_sequential_picker() {
        let group_generators = vec![2.0, 0.5];
        let picker = SequentialPicker::new(&group_generators, 5);

        let mut picker_iter = picker.into_iter();

        assert_eq!(picker_iter.next(), Some(2.0));
        assert_eq!(picker_iter.next(), Some(0.5));
        assert_eq!(picker_iter.next(), Some(2.0));
        assert_eq!(picker_iter.next(), Some(0.5));
        assert_eq!(picker_iter.next(), Some(2.0));
        assert_eq!(picker_iter.next(), None);

        let group_generators = vec![2.0, 0.5, -4.0, -0.25];
        let picker = SequentialPicker::new(&group_generators, 5);

        let mut picker_iter = picker.into_iter();

        assert_eq!(picker_iter.next(), Some(2.0));
        assert_eq!(picker_iter.next(), Some(0.5));
        assert_eq!(picker_iter.next(), Some(-4.0));
        assert_eq!(picker_iter.next(), Some(-0.25));
        assert_eq!(picker_iter.next(), Some(2.0));
        assert_eq!(picker_iter.next(), None);
    }
}
