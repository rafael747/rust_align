//extern crate smith_waterman;
extern crate lib_align;

fn main() {

    let s1 = "CGATTAT";
    let s2 = "GACGTTA";
    let mut smitty = lib_align::SmithWaterman::new(s1.to_string(), s2.to_string());
    let mut needle = lib_align::NeedlemanWunsch::new(s1.to_string(), s2.to_string());

    println!("\nSequencia 1 {} ", s1);
    println!("Sequencia 2 {} \n", s2);

    let local_alignment = smitty.align();
    let global_alignment = needle.align();

    println!("Local alingnment");	
    println!("{:?}", smitty);
    println!("{:?}", local_alignment);

    println!("\nGlobal alingnment");	
    println!("{:?}", needle);
    println!("{:?}", global_alignment);

}



