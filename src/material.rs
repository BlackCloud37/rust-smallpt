#[derive(Debug, Clone, Copy)]
pub enum ReflT {
    DIFF,
    SPEC,
    REFR,
}

pub trait Material {}
