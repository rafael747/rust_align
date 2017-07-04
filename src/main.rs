extern crate lib_align;
use std::io::{self};

fn main() {

//    let s1 = "GAATTCAGTTA";
//    let s2 = "GGATCGA";

    // Leitura da entrada padrão
    let mut s1 = String::new();
    let _ = io::stdin().read_line(&mut s1);

    let mut s2 = String::new();
    let _ = io::stdin().read_line(&mut s2);

    // Remove trailing "\n"
    s1.pop(); 
    s2.pop(); 
  
    let mut smitty = lib_align::SmithWaterman::new(s1.to_string(), s2.to_string());
    let mut needle = lib_align::NeedlemanWunsch::new(s1.to_string(), s2.to_string());

 //   println!("\nSequencia 1 {} ", s1);
 //   println!("Sequencia 2 {} \n", s2);

    let local_alignment = smitty.align();
    let global_alignment = needle.align();

    println!("Local alingnment score");	
 //   println!("{:?}", smitty);
    println!("{:?}", local_alignment.2);

    println!("\nGlobal alingnment score");	
 //   println!("{:?}", needle);
    println!("{:?}", global_alignment.2);

}



