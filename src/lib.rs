extern crate nalgebra;
use nalgebra::DMat;
use std::fmt::{Debug, Formatter, Result};

static MATCH: isize = 5;
static MISMATCH: isize = -3;
static GAP: isize = -4;

/// ESTRUTURA DO OBJETO SMITH-WATERMAN
pub struct SmithWaterman {
    pub sequencia1: Vec<char>,  //Seq1
    pub sequencia2:  Vec<char>,  //Seq2
    pub matrix:  DMat<isize>,
    pub matched: isize,
    pub missed: isize,
    pub gap: isize,
}

// Enumeração dos movimentos no backtrack
pub enum GraphMovements {Blank, Left, Top, Diagonal}

// Implentações das funções do objeto    
impl SmithWaterman{


    // Construtor
    pub fn new(sequencia1: String, sequencia2: String) -> SmithWaterman {
        let matrix = DMat::new_zeros(0,0);
        SmithWaterman{matrix: matrix, sequencia1: sequencia1.chars().collect(),
        sequencia2: sequencia2.chars().collect(), matched: MATCH, missed: MISMATCH, gap: GAP}
    }

    // Função faz o alinhamento e retorna os valores   
    pub fn align(&mut self) -> (String, String, isize){
        let max_point = self.set_matrix_loops();
        return self.local_alignment(max_point);
    }


    // função que calcula novos valores 
    fn penalty(&self, value: isize, penalty_value: isize) -> isize{
        match value.checked_add(penalty_value){
            Some(i) =>{
                if i<0 { 0 }else { i }
            },
            _ => {0}
        }
    }

    //Função que retorna o maior entre os valores possiveis
    fn calculated_movement(&self, row: usize, col: usize) -> isize{
        let left = if col>=1 {self.penalty(self.matrix[(row, col-1)], self.gap)} else {0};
        let top = if row>=1 {self.penalty(self.matrix[(row-1, col)], self.gap)} else {0};
        let diagonal_value = if row>=1 && col>=1 {self.matrix[(row-1, col-1)]} else {0};
        let diagonal_match = if row == 0 || col == 0{
            0
        }else if self.sequencia2.get(row-1).unwrap() == self.sequencia1.get(col-1).unwrap(){
            self.penalty(diagonal_value, self.matched)
        }else{
            self.penalty(diagonal_value, self.missed)
        };
        std::cmp::max(left, std::cmp::max(top, diagonal_match))
    }

    // Função que preenche a matriz e retorna o valor maximo
    pub fn set_matrix_loops(&mut self) -> (usize, usize) {
        self.matrix = nalgebra::DMat::new_zeros(self.sequencia2.len()+1, self.sequencia1.len()+1);
        let mut max_point = (0,0);
        let mut max = 0;
        for row in (1..self.sequencia2.len()+1){
            for col in (1..self.sequencia1.len()+1){
                let n = self.calculated_movement(row, col);
                if n >= max{
                    max = n;
                    max_point = (row,col);
                }
                self.matrix[(row, col)] = n;
            }
        }
        return max_point;

    }


    // Função que faz o backtrack
    fn local_alignment(&self, max_point_value: (usize, usize)) -> (String, String, isize){
        let mut max_point = max_point_value;
        let mut max = self.matrix[max_point];
	    let score = max;
        let mut last_movement = GraphMovements::Blank;
        let mut sequencia1_alignment: Vec<char> =  Vec::with_capacity(self.sequencia1.len());
        let mut sequencia2_alignment: Vec<char> =  Vec::with_capacity(self.sequencia2.len());
        

	    while max > 0 {
            let (row, col) = max_point;
            let one = self.sequencia1.get(col-1).unwrap().clone();
            let two = self.sequencia2.get(row-1).unwrap().clone();
            let top = self.matrix[(row-1, col)];
            let left = self.matrix[(row, col-1)];
            let diagonal = self.matrix[(row-1, col-1)];

	    let diagonal_match = if self.sequencia2.get(row-1).unwrap() == self.sequencia1.get(col-1).unwrap(){
                   self.penalty(diagonal, self.matched)
                }else{
                   self.penalty(diagonal, self.missed)};


	    max_point = if diagonal_match == self.matrix[(row, col)]{
                   last_movement = GraphMovements::Diagonal;
                   (row-1, col-1)
                } else if left == self.matrix[(row, col)]{
                   last_movement = GraphMovements::Left;
                   (row, col-1)
                } else {
                   last_movement = GraphMovements::Top;
                   (row-1, col)
               };
	    
        max = self.matrix[max_point];


        match last_movement {
                GraphMovements::Blank  => {
                    sequencia1_alignment.push(one);
                    sequencia2_alignment.push(two);
                },
                GraphMovements::Diagonal => {
                    sequencia1_alignment.push(one);
                    sequencia2_alignment.push(two);
                },
                GraphMovements::Top => {
                    sequencia1_alignment.push('-');
                    sequencia2_alignment.push(two);
                },
                GraphMovements::Left => {
                    sequencia1_alignment.push(one);
                    sequencia2_alignment.push('-');
                },
            }
        };
        sequencia1_alignment.reverse();
        let x1: String = sequencia1_alignment.into_iter().collect();
        sequencia2_alignment.reverse();
        let x2: String = sequencia2_alignment.into_iter().collect();
        return (x1,x2, score)
    }
}


// Função de Debug para imprimir a matriz

impl Debug for SmithWaterman {
    fn fmt(&self, form:&mut Formatter) -> Result {
        //nrows already has an extra over the sequence counts for the row of zeros
        for row in 0..self.matrix.nrows()+1 {
            for col in 0..self.matrix.ncols()+1 {
                let _ = if col==0 && row>1{
                    write!(form, "{:>5}", self.sequencia2.get(row-2).unwrap().to_string())
                } else if row==0 && col>1{
                    write!(form, "{:>5}", self.sequencia1.get(col-2).unwrap().to_string())
                } else if row>=1 && col>=1{
                    write!(form, "{:>5}", self.matrix[(row-1,col-1)])
                }else{
                    write!(form, "{:>5}", "-")
                };
            }
            let _ = write!(form, "\n");
        }
        write!(form, "\n")
    }
}

///////////////////////////////////////////////////////////////////////////////////////////////


// ESTRUTURA DO OBJETO NEEDLEMAN-WUNSCH
pub struct NeedlemanWunsch {
    pub sequencia1: Vec<char>,
    pub  sequencia2:  Vec<char>,
    pub matrix:  DMat<isize>,
    pub matched: isize,
    pub missed: isize,
    pub gap: isize,
}

// Implementação das funções do objeto
impl NeedlemanWunsch{
   
    // Contrutor
    pub fn new(sequencia1: String, sequencia2: String) -> NeedlemanWunsch {
        let matrix = DMat::new_zeros(0,0);
        NeedlemanWunsch{matrix: matrix, sequencia1: sequencia1.chars().collect(),
        sequencia2: sequencia2.chars().collect(), matched: MATCH, missed: MISMATCH, gap: GAP}
    }

    pub fn align(&mut self) -> (String, String, isize){
        let last_point = self.set_matrix_loops();
        return self.global_alignment(last_point);
    }

    //Função que calcula novos valores
    fn penalty(&self, value: isize, penalty_value: isize) -> isize{
        match value.checked_add(penalty_value){
            Some(i) =>{ i },
            _ => {0}
        }
    }


    fn calculated_movement(&self, row: usize, col: usize) -> isize{
	
        let left = if col>=1 {self.penalty(self.matrix[(row, col-1)], self.gap)} else {0};
        let top = if row>=1 {self.penalty(self.matrix[(row-1, col)], self.gap)} else {0}; 
        let diagonal_value = if row>=1 && col>=1 {self.matrix[(row-1, col-1)]} else {0};
        let diagonal_match = if row == 0 || col == 0{
            0
        }else if self.sequencia2.get(row-1).unwrap() == self.sequencia1.get(col-1).unwrap(){
            self.penalty(diagonal_value, self.matched)
        }else{
            self.penalty(diagonal_value, self.missed) 
        };

        std::cmp::max(left, std::cmp::max(top, diagonal_match))
    }

    
    pub fn set_matrix_loops(&mut self) -> (usize, usize) {
        self.matrix = nalgebra::DMat::new_zeros(self.sequencia2.len()+1, self.sequencia1.len()+1);

	    for i in 0..self.sequencia2.len()+1 {
		    self.matrix[(i,0)]=self.gap * i as isize;
	    }

	    for i in 0..self.sequencia1.len()+1 {
		    self.matrix[(0,i)]=self.gap * i as isize;
	    }

        let last_point = (self.sequencia2.len(), self.sequencia1.len());
        //let mut max = 0;
        for row in (1..self.sequencia2.len()+1){
            for col in (1..self.sequencia1.len()+1){
                let n = self.calculated_movement(row, col);
                self.matrix[(row, col)] = n;
            }
        }
        return last_point;
    }

    fn global_alignment(&self, last_point_value: (usize, usize)) -> (String, String, isize){
        let mut last_point = last_point_value;
        let score = self.matrix[last_point];
        let mut last_movement = GraphMovements::Blank;
        let mut sequencia1_alignment: Vec<char> =  Vec::with_capacity(self.sequencia1.len());
        let mut sequencia2_alignment: Vec<char> =  Vec::with_capacity(self.sequencia2.len());
        

    	while last_point.0 > 0 || last_point.1 > 0 {
            let (row, col) = last_point;
            let one = if col>0 {self.sequencia1.get(col-1).unwrap().clone()} else {'-'};
            let two = if row>0 {self.sequencia2.get(row-1).unwrap().clone()} else {'-'};
          //  let top = if row>0 {self.penalty(self.matrix[(row-1, col)], self.gap)} else {0};  ///unused because in "else"
            let left = if col>0 {self.penalty(self.matrix[(row, col-1)], self.gap)} else {0};
            let diagonal = if row>0 && col>0 {self.matrix[(row-1, col-1)]} else {0};
		
	    if row>0 && col>0{

		let diagonal_match = if self.sequencia2.get(row-1).unwrap() == self.sequencia1.get(col-1).unwrap(){
           	   self.penalty(diagonal, self.matched)
		}else{
            	   self.penalty(diagonal, self.missed)};

	//	last  = std::cmp::max(left, std::cmp::max(top, diagonal_match));

                last_point = if diagonal_match == self.matrix[(row, col)]{
                   last_movement = GraphMovements::Diagonal;
                   (row-1, col-1)
                } else if left == self.matrix[(row, col)]{
                   last_movement = GraphMovements::Left;
                   (row, col-1)
                } else {
                   last_movement = GraphMovements::Top;
                   (row-1, col)
               };

	    } else if row>0{
	//	last  = top;
		last_movement = GraphMovements::Top;
                last_point = (row-1, col);

	    } else{
	//	last  = left;
          	last_movement = GraphMovements::Left;
                last_point = (row, col-1);

	    }

            match last_movement {
                GraphMovements::Blank  => {
                    sequencia1_alignment.push(one);
                    sequencia2_alignment.push(two);
                },
                GraphMovements::Diagonal => {
                    sequencia1_alignment.push(one);
                    sequencia2_alignment.push(two);
                },
                GraphMovements::Top => {
                    sequencia1_alignment.push('-');
                    sequencia2_alignment.push(two);
                },
                GraphMovements::Left => {
                    sequencia1_alignment.push(one);
                    sequencia2_alignment.push('-');
                },
            }
        };
        sequencia1_alignment.reverse();
        let x1: String = sequencia1_alignment.into_iter().collect();
        sequencia2_alignment.reverse();
        let x2: String = sequencia2_alignment.into_iter().collect();

        return (x1,x2,score)
    }
}

impl Debug for NeedlemanWunsch {
    fn fmt(&self, form:&mut Formatter) -> Result {
        //nrows already has an extra over the sequence counts for the row of zeros
        for row in 0..self.matrix.nrows()+1 {
            for col in 0..self.matrix.ncols()+1 {
                let _ = if col==0 && row>1{
                    write!(form, "{:>5}", self.sequencia2.get(row-2).unwrap().to_string())
                } else if row==0 && col>1{
                    write!(form, "{:>5}", self.sequencia1.get(col-2).unwrap().to_string())
                } else if row>=1 && col>=1{
                    write!(form, "{:>5}", self.matrix[(row-1,col-1)])
                }else{
                    write!(form, "{:>5}", "-")
                };
            }
            let _ = write!(form, "\n");
        }
        write!(form, "\n")
    }
}
