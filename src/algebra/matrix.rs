use std::cmp::Ordering;
// use super::group::Semigroup;
use std::ops::{Add, AddAssign, Div, Mul, Neg, Sub};
use std::fmt::Display;

use super::group::{Monoid, Semigroup};

#[derive(Clone)]
pub struct Matrix<T> {
    data: Vec<T>,
    rows: usize,
    cols: usize,
}

impl<T> Display for Matrix<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        for i in 0..self.rows {
            for j in 0..self.cols {
                write!(f, "{} ", self.data[i * self.cols + j])?;
            }
            writeln!(f)?;
        }
        Ok(())
    }
}

impl<T> Matrix<T> {
    pub fn new(data: Vec<T>, rows: usize, cols: usize) -> Self {
        if data.len() != rows * cols {
            panic!("Data length does not match matrix dimensions");
        }
        Matrix { data, rows, cols }
    }
}

impl<T> Mul for &Matrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + std::ops::AddAssign + Copy + Default,
{
    type Output = Matrix<T>;

    fn mul(self, other: Self) -> Matrix<T> {
        if self.cols != other.rows {
            panic!("Matrix dimensions do not match for multiplication");
        }
        let mut result: Matrix<T> = Matrix {
            data: vec![T::default(); self.rows * other.cols],
            rows: self.rows,
            cols: other.cols,
        };
        for i in 0..self.rows {
            for j in 0..other.cols {
                for k in 0..self.cols {
                    result.data[i * other.cols + j] +=
                        self.data[i * self.cols + k] * other.data[k * other.cols + j];
                }
            }
        }
        result
    }
}

impl<T> Add for &Matrix<T>
where
    T: std::ops::Add<Output = T> + Copy + Default,
{
    type Output = Matrix<T>;

    fn add(self, other: Self) -> Self::Output {
        if self.rows != other.rows || self.cols != other.cols {
            panic!("Matrix dimensions do not match for addition");
        }
        let mut result: Matrix<T> = Matrix {
            data: vec![T::default(); self.rows * self.cols],
            rows: self.rows,
            cols: self.cols,
        };
        for i in 0..self.rows {
            for j in 0..self.cols {
                result.data[i * self.cols + j] =
                    self.data[i * self.cols + j] + other.data[i * other.cols + j];
            }
        }
        result
    }
}

#[derive(Clone)]
pub struct SquareMatrix<T> {
    matrix: Matrix<T>,
}

impl<T> SquareMatrix<T> {
    pub fn new(data: Vec<T>, n: usize) -> Self {
        let matrix = Matrix::new(data, n, n);
        SquareMatrix { matrix }
    }

    pub fn size(&self) -> usize {
        self.matrix.rows
    }

    pub fn into_matrix(self) -> Matrix<T> {
        self.matrix
    }
}

impl<T> TryFrom<Matrix<T>> for SquareMatrix<T> {
    type Error = String;

    fn try_from(matrix: Matrix<T>) -> Result<Self, Self::Error> {
        if matrix.rows == matrix.cols {
            Ok(SquareMatrix { matrix })
        } else {
            Err("Matrix is not square".to_string())
        }
    }
}

impl<T> Mul for &SquareMatrix<T>
where
    T: Mul<Output = T> + Add<Output = T> + AddAssign + Copy + Default,
{
    type Output = SquareMatrix<T>;

    fn mul(self, other: Self) -> Self::Output {
        let product_matrix = &self.matrix * &other.matrix;
        SquareMatrix {
            matrix: product_matrix,
        }
    }
}

impl<T> Add for &SquareMatrix<T>
where
    T: Add<Output = T> + Copy + Default,
{
    type Output = SquareMatrix<T>;

    fn add(self, other: Self) -> Self::Output {
        let sum_matrix = &self.matrix + &other.matrix;
        SquareMatrix {
            matrix: sum_matrix,
        }
    }
}

impl<T> Display for SquareMatrix<T>
where
    T: Display,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.matrix.fmt(f)
    }
}

pub trait Zero {
    fn zero() -> Self;
}

pub trait One {
    fn one() -> Self;
}

// 为常见类型实现特征

// Define the macro `impl_num_traits` if it is missing
macro_rules! impl_num_traits {
    ($($t:ty),*) => {
        $(
            impl Zero for $t {
                fn zero() -> Self {
                    0 as $t
                }
            }

            impl One for $t {
                fn one() -> Self {
                    1 as $t
                }
            }
        )*
    };
}

impl_num_traits!(i32, i64, f32, f64);


// 方阵实现
impl<T> SquareMatrix<T> {
    /// 计算方阵的迹（对角线元素之和）
    pub fn trace(&self) -> T 
    where
        T: Zero + Add<Output = T> + Clone,
    {
        let n = self.size();
        let mut trace = T::zero();
        for i in 0..n {
            trace = trace + self.matrix.data[i * n + i].clone();
        }
        trace
    }

    /// 生成单位矩阵
    pub fn identity(n: usize) -> Self 
    where
        T: One + Zero + Clone,
    {
        let mut data = Vec::with_capacity(n * n);
        for i in 0..n {
            for j in 0..n {
                data.push(if i == j { T::one() } else { T::zero() });
            }
        }
        SquareMatrix::new(data, n)
    }

    /// 递归计算行列式
    pub fn determinant(&self) -> T 
    where
        T: Clone + Zero + One + Add<Output = T> + Sub<Output = T> + Mul<Output = T> + Neg<Output = T>,
    {
        let n = self.size();
        match n {
            1 => self.matrix.data[0].clone(),
            2 => {
                let a = &self.matrix.data[0];
                let b = &self.matrix.data[1];
                let c = &self.matrix.data[2];
                let d = &self.matrix.data[3];
                a.clone() * d.clone() - b.clone() * c.clone()
            }
            _ => {
                let mut det = T::zero();
                for col in 0..n {
                    let sign = if col % 2 == 0 { T::one() } else { -T::one() };
                    det = det + sign * self.matrix.data[col].clone() * self.minor(0, col).determinant();
                }
                det
            }
        }
    }

    /// 高斯-约旦消元法求逆矩阵
    pub fn inverse(&self) -> Result<Self, String> 
    where
        T: Clone + Default + PartialEq + PartialOrd 
            + Add<Output = T> + Sub<Output = T> 
            + Mul<Output = T> + Div<Output = T> + Neg<Output = T> 
            + Zero + One,
    {
        let n = self.size();
        let mut mat = self.matrix.clone();
        let mut aug: Vec<T> = SquareMatrix::identity(n).matrix.data;

        for i in 0..n {
            // 寻找主元
            let pivot = (i..n).max_by(|&a, &b| 
                mat.data[a * n + i].partial_cmp(&mat.data[b * n + i])
                    .unwrap_or(Ordering::Equal)
                ).unwrap();

            if mat.data[pivot * n + i] == T::zero() {
                return Err("Matrix is singular".into());
            }

            // 交换行
            if pivot != i {
                for j in 0..n {
                    mat.data.swap(i * n + j, pivot * n + j);
                    aug.swap(i * n + j, pivot * n + j);
                }
            }

            // 归一化当前行
            let diag = mat.data[i * n + i].clone();
            for j in 0..n {
                mat.data[i * n + j] = mat.data[i * n + j].clone() / diag.clone();
                aug[i * n + j] = aug[i * n + j].clone() / diag.clone();
            }

            // 消去其他行
            for k in 0..n {
                if k != i {
                    let factor = mat.data[k * n + i].clone();
                    for j in 0..n {
                        mat.data[k * n + j] = mat.data[k * n + j].clone() - factor.clone() * mat.data[i * n + j].clone();
                        aug[k * n + j] = aug[k * n + j].clone() - factor.clone() * aug[i * n + j].clone();
                    }
                }
            }
        }

        Ok(SquareMatrix::new(aug, n))
    }

    /// 获取子矩阵
    fn minor(&self, row: usize, col: usize) -> Self 
    where
        T: Clone,
    {
        let n = self.size();
        let mut data = Vec::with_capacity((n-1)*(n-1));
        for i in 0..n {
            for j in 0..n {
                if i != row && j != col {
                    data.push(self.matrix.data[i * n + j].clone());
                }
            }
        }
        SquareMatrix::new(data, n-1)
    }
}