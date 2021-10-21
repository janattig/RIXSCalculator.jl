# TEST FOR THE PHYSICS BEHIND THE CODE
# it mainly tests eigenenergies with specific hamiltonians and consistency in various instances


@testset "Physics tests" begin
    
    @testset "Sample system tests" begin
    
        @testset "BdotS only" begin
            
            @testset "1h, 1s" begin
                
                @testset "1h, 1s: single-particle picture" begin

                    h=1
                    s=1
                    Bstr=1.0
                    Bdir=generate_rvos()
                    B=Bstr*Bdir

                    # basis construction
                    basis_sp=getT2GBasisLS()        

                    # hamiltonian construction
                    hamiltonian=MagneticFieldOperator(basis_sp, Bstr, Bdir)

                    # obtain eigensystem
                    es=eigensystem(hamiltonian)

                    @test size(matrix_representation(hamiltonian))==(6,6)
                    @test length(es[:values])==6
                    @test length(es[:vectors])==6

                    # correct results:
                    Epm=norm(B)/2
                    E1=-Epm
                    E2=Epm

                    #test
                    @test abs.(es[:values]-[E1,E1,E1,E2,E2,E2])<1e-6*ones(6)

                end

                @testset "1h, 1s: multi-particle picture" begin

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

                    @test size(matrix_representation(hamiltonian))==(6,6)
                    @test length(es[:values])==6
                    @test length(es[:vectors])==6

                    # correct results:
                    Epm=norm(B)/2
                    E1=-Epm
                    E2=Epm

                    #test
                    @test abs.(es[:values]-[E1,E1,E1,E2,E2,E2])<1e-6*ones(6)

                end
                
            end # end 1h, 1s
            
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
                
                
                # consistency tests
                @test size(matrix_representation(hamiltonian))==(15,15)
                @test length(es[:values])==15
                @test length(es[:vectors])==15
                
                
                # correct results:
                E1=-Bstr
                E2=0.0
                E3=Bstr
                energies=zeros(15)
                energies[1:3]=E1*ones(3)
                energies[13:15]=E3*ones(3)

                # test 
                @test abs.(es[:values]-energies)<1e-6*ones(15)
                
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
                
                # consistency tests
                @test size(matrix_representation(hamiltonian))==(20,20)
                @test length(es[:values])==20
                @test length(es[:vectors])==20

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
                @test abs.(es[:values]-energies)<1e-6*ones(20)

            end
            
            @testset "1h, 2s, BdotS on one site only" begin
                
                h=1
                s=2
                Bstr=1.0
                Bdir=generate_rvos()
                B=Bstr*Bdir

                # basis construction
                basis_sp=getT2GBasisLS()        
                basis_ms=getMultiSiteBasis(basis_sp,s)

                # hamiltonian construction
                hamiltonian=MagneticFieldOperator(basis_ms, 1, Bstr, Bdir)

                # obtain eigensystem
                es=eigensystem(hamiltonian)

                @test size(matrix_representation(hamiltonian))==(12,12)
                @test length(es[:values])==12
                @test length(es[:vectors])==12

                # correct results:
                energies=zeros(12)
                energies[1:3]=-(Bstr/2)*ones(3)
                energies[4:9]=zeros(6)
                energies[10:12]=(Bstr/2)*ones(3)

                # #test
                @test abs.(es[:values]-energies)<1e-6*ones(12)
                
            end
            
            
            end # end BdotS only
        
        @testset "LdotS only" begin
            
            @testset "1h, 1s" begin
                
                @testset "1h, 1s: single-particle picture" begin

                    h=1
                    s=1
                    lambda=3.0

                    # basis construction
                    basis_sp=getT2GBasisLS()        

                    # hamiltonian construction
                    hamiltonian=SpinOrbitOperator(basis_sp, lambda)

                    # obtain eigensystem
                    es=eigensystem(hamiltonian)
                    es[:values]

                    # consistency tests
                    @test size(matrix_representation(hamiltonian))==(6,6)
                    @test length(es[:values])==6
                    @test length(es[:vectors])==6

                    # correct results:
                    E1=-lambda
                    E2=+lambda/2
                    energies=zeros(length(es[:values]))
                    energies[1:2]=E1*ones(2)
                    energies[3:6]=E2*ones(4)

                    #test
                    @test abs.(es[:values]-energies)<1e-6*ones(length(es[:values]))

                end

                @testset "1h, 1s: multi-particle picture" begin

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

                    # consistency tests
                    @test size(matrix_representation(hamiltonian))==(6,6)
                    @test length(es[:values])==6
                    @test length(es[:vectors])==6

                    # correct results:
                    E1=-lambda
                    E2=+lambda/2
                    energies=zeros(length(es[:values]))
                    energies[1:2]=E1*ones(2)
                    energies[3:6]=E2*ones(4)

                    #test
                    @test abs.(es[:values]-energies)<1e-6*ones(length(es[:values]))

                end
            end # end 1h, 1s 
            
            @testset "2h, 1s" begin
                
                h=2
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
                
                # consistency tests
                @test size(matrix_representation(hamiltonian))==(15,15)
                @test length(es[:values])==15
                @test length(es[:vectors])==15

                # correct results:
                E1=-2*lambda
                E2=-lambda/2
                E3=lambda

                energies=zeros(length(es[:values]))
                energies[1]=E1
                energies[2:8]=E2*ones(7)
                energies[9:15]=E3*ones(7)

                #test
                @test abs.(es[:values]-energies)<1e-6*ones(length(es[:values])) 
                
            end
            
            @testset "3h, 1s" begin
                
                h=3
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

                # consistency tests
                @test size(matrix_representation(hamiltonian))==(20,20)
                @test length(es[:values])==20
                @test length(es[:vectors])==20

                # correct results:
                E1=-3*lambda/2
                E2=0.0
                E3= 3*lambda/2

                energies=zeros(length(es[:values]))
                energies[1:4]=E1*ones(4)
                energies[5:16]=E2*ones(12)
                energies[17:20]=E3*ones(4)

                #test
                @test abs.(es[:values]-energies)<1e-6*ones(length(es[:values])) 
                
            end
            
            @testset "1h, 2s: LdotS on one site only" begin
                
                h=1
                s=2
                lambda=3.0

                # basis construction
                basis_sp=getT2GBasisLS()        
                basis_ms=getMultiSiteBasis(basis_sp,s)


                # hamiltonian construction
                hamiltonian=SpinOrbitOperator(basis_ms, 1, lambda) #+ SpinOrbitOperator(basis_ms, 2, lambda)

                # obtain eigensystem
                es=eigensystem(hamiltonian)

                # consistency tests
                @test size(matrix_representation(hamiltonian))==(12,12)
                @test length(es[:values])==12
                @test length(es[:vectors])==12

                # correct results:
                E1=-lambda
                E2=+lambda/2
                energies=zeros(12)
                energies[1:2]=E1*ones(2)
                energies[3:8]=zeros(6)
                energies[9:12]=E2*ones(4)

                #test
                @test abs.(es[:values]-energies)<1e-6*ones(length(es[:values]))
                
            end
            
        end
    
    end

    
    @testset "Literature tests" begin
        
        @testset "Comparison with Ament et al." begin #Khaliullin and van der Brink
            
            # auxiliary functions
            function theta(lambda, Delta)
                return 0.5*atan( 2*sqrt(2)*lambda/(lambda-2*Delta) )
            end
            function ament_energies(lambda, Delta, theta)
                Ef=lambda/( sqrt(2)*tan(theta) )
                Eg=-Delta-lambda/2
                Eh=-( lambda*tan(theta) )/sqrt(2)
                # make sure to return sort(list)
                return sort([Ef,Eh,Eg])
            end
            
            # choice of parameters
            lambda=rand(-10.0:10.0)
            Delta=rand(-10.0:10.0)

            # basis construction
            basis_sp=getT2GBasisXYZ()
            basis_mp=getMultiParticleBasis(getMultiSiteBasis(basis_sp,1),1)

            # hamiltonian construction
            hamiltonian=DistortionOperator(basis_mp,1,-Delta, [0,0,1]) + SpinOrbitOperator(basis_mp, 1, -lambda)
            es=eigensystem(hamiltonian)

            # ament solution
            energies=zeros(6)
            E1,E2,E3=ament_energies(lambda, Delta, theta(lambda,Delta))
            energies[1:2]=E1*ones(2)
            energies[3:4]=E2*ones(2)
            energies[5:6]=E3*ones(2)

            #test
            @test abs.(es[:values]-energies)<1e-6*ones(6)
            
        end
        
        
    end
    
    
    
end #end physics tests