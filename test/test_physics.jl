# TEST FOR THE PHYSICS BEHIND THE CODE
# it mainly tests eigenenergies with specific hamiltonians and consistency in various instances

# Auxiliary functions:
"""
    generate_rvos()

Generates a Random Vector On a Sphere.
"""
function generate_rvos()
    theta=rand()*2*pi    #rand(Uniform(0.0,2*pi))
    z=-1.0+2*rand()    #rand(Uniform(0.0,1.0))
    
    x=sqrt(1-z^2)*cos(theta)
    y=sqrt(1-z^2)*sin(theta)
    
    return [x,y,z]
end



@testset "Physics tests" begin
    
    @testset "Eigenvalue tests" begin
    
        @testset "BdotS only" begin
            
            @testset "1h, 1s" begin
               
                h=1
                s=1
                Bstr=1.0
                Bdir=generate_rvos()
                B=Bstr*Bdir
                
                # basis construction
                basis_sp=getT2GBasisLS()        
                basis_ms=getMultiSiteBasis(basis_sp,s)
                basis_mp=getMultiParticleBasis(basis_ms,h)
                
                # hamiltonian construction
                hamiltonian=MagneticFieldOperator(basis_mp, 1, Bstr, Bdir)
                
                # obtain eigensystem
                es=eigensystem(hamiltonian)
                
                # correct results:
                Epm=norm(B)/2
                E1=-Epm
                E2=Epm
                
                #test
                @test (es[:values]-[E1,E1,E1,E2,E2,E2])<1e-6*ones(6)
                
            end
            
            @testset "2h, 1s" begin
               
                h=2
                s=1
                Bstr=1.0
                Bdir=generate_rvos()
                B=Bstr*Bdir

                # basis construction
                basis_sp=getT2GBasisLS()        
                basis_ms=getMultiSiteBasis(basis_sp,s)
                basis_mp=getMultiParticleBasis(basis_ms,h)

                # hamiltonian construction
                hamiltonian=MagneticFieldOperator(basis_mp, 1, Bstr, Bdir)

                # correct results:
                E1=-1.0
                E2=0.0
                E3=1.0
                energies=zeros(15)
                energies[1:3]=E1*ones(3)
                energies[13:15]=E3*ones(3)

                # obtain eigensystem
                es=eigensystem(hamiltonian)
                es[:values]

                # test 
                @test (es[:values]-energies)<1e-6*ones(15)
                
            end
            
            
        end
    
    end


    @testset "Eigenvector tests" begin
        
        
        
    end
    
    @testset "Literature tests" begin
        
        
        
    end
    
    
    
end #end physics tests