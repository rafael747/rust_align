//extern crate smith_waterman;
extern crate lib_align;

fn main() {

    let s1 = "CGATTAT";
    let s2 = "GACGTTA";
    let mut smitty = lib_align::SmithWaterman::new(s1.to_string(), s2.to_string());

    println!("\nsequencia 1 {} ", s1);
    println!("sequencia 2 {} \n", s2);

    let alignment = smitty.align();
    println!("{:?}", smitty);
    println!("{:?}", alignment);

}



