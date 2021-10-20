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

                # obtain eigensystem
                es=eigensystem(hamiltonian)
                es[:values]
                
                
                # correct results:
                E1=-Bstr
                E2=0.0
                E3=Bstr
                energies=zeros(15)
                energies[1:3]=E1*ones(3)
                energies[13:15]=E3*ones(3)

                # test 
                @test (es[:values]-energies)<1e-6*ones(15)
                
            end
            
            @testset "3h, 1s" begin
               
                h=3
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
                es[:values]

                # correct results:
                E1=-3*Bstr/2
                E2=-Bstr/2
                E3= Bstr/2
                E4=3*Bstr/2

                energies=zeros(20)
                energies[1]=E1
                energies[2:10]=E2*ones(9)
                energies[11:19]=E3*ones(9)
                energies[20]=E4

                # test 
                @test (es[:values]-energies)<1e-6*ones(20)

            end
            
            
            end # end BdotS only
        
        @testset "LdotS only" begin
            
            @testset "1h, 1s" begin
               
                h=1
                s=1
                lambda=3.0

                # basis construction
                basis_sp=getT2GBasisLS()        
                basis_ms=getMultiSiteBasis(basis_sp,s)
                basis_mp=getMultiParticleBasis(basis_ms,h)

                # hamiltonian construction
                hamiltonian=SpinOrbitOperator(basis_mp, 1, lambda)

                # obtain eigensystem
                es=eigensystem(hamiltonian)
                es[:values]

                # correct results:
                E1=-lambda
                E2=+lambda/2
                energies=zeros(length(es[:values]))
                energies[1:2]=E1*ones(2)
                energies[3:6]=E2*ones(4)

                #test
                @test (es[:values]-energies)<1e-6*ones(length(es[:values]))
                
            end
            
        end
    
    end


    @testset "Eigenvector tests" begin
        
        
        
    end
    
    @testset "Literature tests" begin
        
        
        
    end
    
    
    
end #end physics tests