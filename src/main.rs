mod algebra;

use algebra::matrix::{Matrix, SquareMatrix};


fn main() -> Result<(), String> {
    let a: Matrix<i32> = Matrix::new(vec![1, 2, 3, 4], 2, 2);
    let b = Matrix::new(vec![5, 6, 7, 8], 2, 2);
    let c = &a * &b;
    let d = &a + &b;
    let e = SquareMatrix::try_from(d.clone())?;
    let f = e.trace();
    println!("Result: {}", c);
    println!("Result: {}", d);
    println!("Result: {}", f);
    Ok(())
}
