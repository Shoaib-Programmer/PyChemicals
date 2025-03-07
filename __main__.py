from .chemicals import Acid, Base

def main():
    # Example usage
    print("Testing PyChemicals module...")
    
    # Create and test an acid
    hcl = Acid(name="Hydrochloric Acid", concentration=0.1)
    print(f"HCl pH: {hcl.pH():.2f}")
    
    # Create and test a base
    naoh = Base(name="Sodium Hydroxide", concentration=0.1)
    print(f"NaOH pH: {naoh.pH():.2f}")

if __name__ == "__main__":
    main()
