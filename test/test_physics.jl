# TEST FOR THE PHYSICS BEHIND THE CODE
# it mainly tests eigenenergies with specific hamiltonians and consistency in various instances


@testset "Physics tests" begin
    
    @testset "Sample system tests" begin
    
        @testset "BdotS only" begin
            
            @testset "1h, 1s" begin
                
                @testset "1h, 1s: single-particle picture" begin

                    h=1
                    s=1
                    particle_type=:hole
                    
                    Bstr=1.0
                    Bdir=generate_rvos()
                    B=Bstr*Bdir

                    # basis construction
                    basis_sp=getT2GBasisLS()        

                    # hamiltonian construction
                    hamiltonian=MagneticFieldOperator(basis_sp, Bstr, Bdir)

                    # obtain eigensystem
                    es=eigensystem(hamiltonian)

                    @test size(matrix_representation(hamiltonian)) == (6,6)
                    @test length(es[:values]) == 6
                    @test length(es[:vectors]) == 6

                    # correct results:
                    Epm=norm(B)/2
                    E1=-Epm
                    E2=Epm

                    #test
                    @test (abs.(es[:values]-[E1,E1,E1,E2,E2,E2]).<1e-6*ones(6)) == trues(6)

                end

                @testset "1h, 1s: multi-particle picture" begin

                    h=1
                    s=1
                    particle_type=:hole
                    
                    Bstr=1.0
                    Bdir=generate_rvos()
                    B=Bstr*Bdir

                    # basis construction
                    basis_sp=getT2GBasisLS()        
                    basis_ms=getMultiSiteBasis(basis_sp,s)
                    basis_mp=getMultiParticleBasis(basis_ms,h)

                    # hamiltonian construction
                    hamiltonian=MagneticFieldOperator(basis_mp, 1, Bstr, Bdir; particle_type=particle_type)

                    # obtain eigensystem
                    es=eigensystem(hamiltonian)

                    @test size(matrix_representation(hamiltonian))==(6,6)
                    @test length(es[:values])==6
                    @test length(es[:vectors])==6

                    # correct results:
                    Epm=norm(B)/2
                    E1=-Epm
                    E2=Epm
                    energies=[E1,E1,E1,E2,E2,E2]

                    #test
                    @test (abs.(es[:values]-[E1,E1,E1,E2,E2,E2]).<1e-6*ones(length(energies))) == trues(length(energies))

                end
                
            end # end 1h, 1s
            
            @testset "1e, 1s" begin

                @testset "1e, 1s: multi-particle picture" begin

                    e=1
                    s=1
                    particle_type=:electron
                    
                    Bstr=1.0
                    Bdir=generate_rvos()
                    B=Bstr*Bdir

                    # basis construction
                    basis_sp=getT2GBasisLS()        
                    basis_ms=getMultiSiteBasis(basis_sp,s)
                    basis_mp=getMultiParticleBasis(basis_ms,e)

                    # hamiltonian construction
                    hamiltonian=MagneticFieldOperator(basis_mp, 1, Bstr, Bdir; particle_type=particle_type)

                    # obtain eigensystem
                    es=eigensystem(hamiltonian)

                    @test size(matrix_representation(hamiltonian))==(6,6)
                    @test length(es[:values])==6
                    @test length(es[:vectors])==6

                    # correct results:
                    Epm=norm(B)/2
                    E1=-Epm
                    E2=Epm
                    energies=[E1,E1,E1,E2,E2,E2]
                    #test
                    @test (abs.(es[:values]-energies).<1e-6*ones(length(es[:values]))) == trues(length(energies))

                end
                
            end # end 1e, 1s
            
            @testset "2h, 1s" begin
               
                h=2
                s=1
                particle_type=:hole
                
                Bstr=1.0
                Bdir=generate_rvos()
                B=Bstr*Bdir

                # basis construction
                basis_sp=getT2GBasisLS()        
                basis_ms=getMultiSiteBasis(basis_sp,s)
                basis_mp=getMultiParticleBasis(basis_ms,h)

                # hamiltonian construction
                hamiltonian=MagneticFieldOperator(basis_mp, 1, Bstr, Bdir; particle_type=particle_type)

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
                @test (abs.(es[:values]-energies).<1e-6*ones(length(es[:values]))) == trues(length(energies))
                
            end
            
            @testset "2e, 1s" begin
               
                e=2
                s=1
                particle_type=:electron
                
                Bstr=1.0
                Bdir=generate_rvos()
                B=Bstr*Bdir

                # basis construction
                basis_sp=getT2GBasisLS()        
                basis_ms=getMultiSiteBasis(basis_sp,s)
                basis_mp=getMultiParticleBasis(basis_ms,e)

                # hamiltonian construction
                hamiltonian=MagneticFieldOperator(basis_mp, 1, Bstr, Bdir; particle_type=particle_type)

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
                @test (abs.(es[:values]-energies).<1e-6*ones(length(es[:values]))) == trues(length(energies))
                
            end
            
            @testset "3h, 1s" begin
               
                h=3
                s=1
                particle_type=:hole
                
                Bstr=1.0
                Bdir=generate_rvos()
                B=Bstr*Bdir

                # basis construction
                basis_sp=getT2GBasisLS()        
                basis_ms=getMultiSiteBasis(basis_sp,s)
                basis_mp=getMultiParticleBasis(basis_ms,h)

                # hamiltonian construction
                hamiltonian=MagneticFieldOperator(basis_mp, 1, Bstr, Bdir; particle_type=particle_type)

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
                @test (abs.(es[:values]-energies).<1e-6*ones(length(es[:values]))) == trues(length(energies))

            end
            
            @testset "3e, 1s" begin
               
                e=3
                s=1
                particle_type=:electron
                
                Bstr=1.0
                Bdir=generate_rvos()
                B=Bstr*Bdir

                # basis construction
                basis_sp=getT2GBasisLS()        
                basis_ms=getMultiSiteBasis(basis_sp,s)
                basis_mp=getMultiParticleBasis(basis_ms,e)

                # hamiltonian construction
                hamiltonian=MagneticFieldOperator(basis_mp, 1, Bstr, Bdir; particle_type=particle_type)

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
                @test (abs.(es[:values]-energies).<1e-6*ones(length(es[:values]))) == trues(length(energies))

            end
            
            @testset "1h, 2s, BdotS on one site only" begin
                
                h=1
                s=2
                particle_type=:hole
                
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
                @test (abs.(es[:values]-energies).<1e-6*ones(length(es[:values]))) == trues(length(energies))
                
            end
            
            
        end # end BdotS only
        
        @testset "LdotS only" begin
            
            @testset "1h, 1s" begin
                
                @testset "1h, 1s: single-particle picture" begin

                    h=1
                    s=1
                    particle_type=:hole
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
                    @test (abs.(es[:values]-energies).<1e-6*ones(length(es[:values]))) == trues(length(energies))

                end

                @testset "1h, 1s: multi-particle picture" begin

                    h=1
                    s=1
                    particle_type=:hole
                    
                    lambda=3.0

                    # basis construction
                    basis_sp=getT2GBasisLS()        
                    basis_ms=getMultiSiteBasis(basis_sp,s)
                    basis_mp=getMultiParticleBasis(basis_ms,h)

                    # hamiltonian construction
                    hamiltonian=SpinOrbitOperator(basis_mp, 1, lambda; particle_type=particle_type)

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
                    @test (abs.(es[:values]-energies).<1e-6*ones(length(es[:values]))) == trues(length(energies))

                end
            end # end 1h, 1s 
            
            @testset "1e, 1s" begin

                @testset "1e, 1s: multi-particle picture" begin

                    e=1
                    s=1
                    particle_type=:electron
                    lambda=3.0

                    # basis construction
                    basis_sp=getT2GBasisLS()        
                    basis_ms=getMultiSiteBasis(basis_sp,s)
                    basis_mp=getMultiParticleBasis(basis_ms,e)

                    # hamiltonian construction
                    hamiltonian=SpinOrbitOperator(basis_mp, 1, lambda; particle_type=particle_type)

                    # obtain eigensystem
                    es=eigensystem(hamiltonian)
                    es[:values]

                    # consistency tests
                    @test size(matrix_representation(hamiltonian))==(6,6)
                    @test length(es[:values])==6
                    @test length(es[:vectors])==6

                    # correct results:
                    E1=-lambda/2
                    E2=+lambda
                    energies=zeros(length(es[:values]))
                    energies[1:4]=E1*ones(4)
                    energies[5:6]=E2*ones(2)

                    #test
                    @test (abs.(es[:values]-energies).<1e-6*ones(length(es[:values]))) == trues(length(energies))

                end
            end # end 1e, 1s 
            
            @testset "2h, 1s" begin
                
                h=2
                s=1
                particle_type=:hole
                lambda=3.0

                # basis construction
                basis_sp=getT2GBasisLS()        
                basis_ms=getMultiSiteBasis(basis_sp,s)
                basis_mp=getMultiParticleBasis(basis_ms,h)

                # hamiltonian construction
                hamiltonian=SpinOrbitOperator(basis_mp, 1, lambda; particle_type=particle_type)

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
                energies[2:9]=E2*ones(8)
                energies[10:15]=E3*ones(6)

                #test
                @test (abs.(es[:values]-energies).<1e-6*ones(length(es[:values]))) == trues(length(energies))
                
            end
            
            @testset "2e, 1s" begin
                
                h=2
                s=1
                particle_type=:electron
                lambda=3.0

                # basis construction
                basis_sp=getT2GBasisLS()        
                basis_ms=getMultiSiteBasis(basis_sp,s)
                basis_mp=getMultiParticleBasis(basis_ms,h)

                # hamiltonian construction
                hamiltonian=SpinOrbitOperator(basis_mp, 1, lambda; particle_type=particle_type)

                # obtain eigensystem
                es=eigensystem(hamiltonian)
                es[:values]
                
                # consistency tests
                @test size(matrix_representation(hamiltonian))==(15,15)
                @test length(es[:values])==15
                @test length(es[:vectors])==15

                # correct results:
                E1=-lambda
                E2=lambda/2
                E3=2*lambda

                energies=zeros(length(es[:values]))
                energies[1:6]=E1*ones(6)
                energies[7:14]=E2*ones(8)
                energies[15]=E3

                #test
                @test (abs.(es[:values]-energies).<1e-6*ones(length(es[:values]))) == trues(length(energies))
                
            end
            
            @testset "3h, 1s" begin
                
                h=3
                s=1
                particle_type=:hole
                lambda=3.0

                # basis construction
                basis_sp=getT2GBasisLS()        
                basis_ms=getMultiSiteBasis(basis_sp,s)
                basis_mp=getMultiParticleBasis(basis_ms,h)

                # hamiltonian construction
                hamiltonian=SpinOrbitOperator(basis_mp, 1, lambda; particle_type=particle_type)

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
                @test (abs.(es[:values]-energies).<1e-6*ones(length(es[:values]))) == trues(length(energies))
                
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
                @test (abs.(es[:values]-energies).<1e-6*ones(length(energies)) ) == trues(length(energies))
                
            end
            
        end # end LdotS tests
        
#         @testset "Hopping" begin
            
        
#             @testset "2s, 1h, multisite picture"
            
            
            
#             end
            
            
#         end #end hopping tests
    
    end


    
    @testset "Literature tests" begin
        
        @testset "Comparison with Ament et al., holes" begin #Khaliullin and van der Brink
            
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
            
            s=1
            p=1
            particle_type=:hole

            # basis construction
            basis_sp=getT2GBasisXYZ()
            basis_mp=getMultiParticleBasis(getMultiSiteBasis(basis_sp,s),p)

            # hamiltonian construction
            hamiltonian=DistortionOperator(basis_mp,1,-Delta, [0,0,1]; particle_type=particle_type) + 
                        SpinOrbitOperator(basis_mp, 1,-lambda; particle_type=particle_type)
            es=eigensystem(hamiltonian)

            # ament solution
            energies=zeros(6)
            E1,E2,E3=ament_energies(lambda, Delta, theta(lambda,Delta))
            energies[1:2]=E1*ones(2)
            energies[3:4]=E2*ones(2)
            energies[5:6]=E3*ones(2)

            #test
            @test (abs.(es[:values]-energies).<1e-6*ones(length(energies)) ) == trues(length(energies))
            
        end
        
        @testset "Comparison with Ament et al., electrons" begin #Khaliullin and van der Brink
            
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
            
            s=1
            p=5
            particle_type=:electron

            # basis construction
            basis_sp=getT2GBasisXYZ()
            basis_mp=getMultiParticleBasis(getMultiSiteBasis(basis_sp,s),p)

            # hamiltonian construction
            hamiltonian=DistortionOperator(basis_mp,1,-Delta, [0,0,1]; particle_type=particle_type) + 
                        SpinOrbitOperator(basis_mp, 1,-lambda; particle_type=particle_type)
            es=eigensystem(hamiltonian)

            # ament solution
            energies=zeros(6)
            E1,E2,E3=ament_energies(lambda, Delta, theta(lambda,Delta))
            energies[1:2]=E1*ones(2)
            energies[3:4]=E2*ones(2)
            energies[5:6]=E3*ones(2)

            #test
            @test (abs.(es[:values]-energies).<1e-6*ones(length(energies)) ) == trues(length(energies))
            
        end
        
        @testset "Comparison with Perkins, Sizyuk and Woelfle" begin
            
            # analytical results
            E1h(U2) = 10*U2
            E0h(U2)=15*U2
            E2h_1(U2,J_H)=6*U2-J_H
            E2h_0(U2,J_H)=6*U2+J_H
            E2h_00(U2,J_H)=6*U2+4*J_H
            
            @testset "0h, 1s" begin
                
                # parameters
                h=0
                s=1

                U1=rand(-1.0:10.0)
                J_H=rand(-1.0:10.0)
                U2=U1-2*J_H

                # basis construction
                basis_sp=getT2GBasisXYZ()
                basis_mp=getMultiParticleBasis(getMultiSiteBasis(basis_sp,s),h)

                # hamiltonian construction
                hamiltonian= MPHolePerkinsWoelfleHamiltonian(basis_mp, 1, U1,U2,J_H)

                # eigensystem
                es=eigensystem(hamiltonian)
                @test length(es[:values])==1
    
                @test abs(es[:values][1]-E0h(U2))<1e-6
                
            end
            
            @testset "1h, 1s" begin
                
                # parameters
                h=1
                s=1

                U1=rand(-1.0:10.0)
                J_H=rand(-1.0:10.0)
                U2=U1-2*J_H

                # basis construction
                basis_sp=getT2GBasisXYZ()
                basis_mp=getMultiParticleBasis(getMultiSiteBasis(basis_sp,s),h)

                # hamiltonian construction
                hamiltonian= MPHolePerkinsWoelfleHamiltonian(basis_mp, 1, U1,U2,J_H)

                # eigensystem
                es=eigensystem(hamiltonian)
                @test length(es[:values])==6
    
                @test (abs.(es[:values]-E1h(U2)*ones(6)).<1e-6*ones(6)) == trues(6)
                
            end
            
            @testset "2h, 1s" begin
                
               # parameters
                h=2
                s=1

                U1=rand(-1.0:10.0)
                J_H=rand(-1.0:10.0)
                U2=U1-2*J_H

                # basis construction
                basis_sp=getT2GBasisXYZ()
                basis_mp=getMultiParticleBasis(getMultiSiteBasis(basis_sp,s),h)

                # hamiltonian construction
                hamiltonian= MPHolePerkinsWoelfleHamiltonian(basis_mp, 1, U1,U2,J_H)

                # eigensystem
                es=eigensystem(hamiltonian) 
                @test length(es[:values])==15
                
                # analytical solutions
                energies=zeros(15)
                energies[1:9]=E2h_1(U2,J_H)*ones(9)
                energies[10:14]=E2h_0(U2,J_H)*ones(5)
                energies[15]=E2h_00(U2,J_H)
                
                # test
                @test (abs.(es[:values]-energies).<1e-6*ones(length(energies)) ) == trues(length(energies)) 
                
            end
            
            
            @testset "6e, 1s" begin
                
                # parameters
                e=6
                s=1

                U1=rand(-1.0:10.0)
                J_H=rand(-1.0:10.0)
                U2=U1-2*J_H

                # basis construction
                basis_sp=getT2GBasisXYZ()
                basis_mp=getMultiParticleBasis(getMultiSiteBasis(basis_sp,s),e)

                # hamiltonian construction
                hamiltonian= MPElectronPerkinsWoelfleHamiltonian(basis_mp, s, U1,U2,J_H)

                # eigensystem
                es=eigensystem(hamiltonian)
                @test length(es[:values])==1
    
                @test abs(es[:values][1]-E0h(U2))<1e-6
                
            end
            
            @testset "5e, 1s" begin
                
                # parameters
                e=5
                s=1

                U1=rand(-1.0:10.0)
                J_H=rand(-1.0:10.0)
                U2=U1-2*J_H

                # basis construction
                basis_sp=getT2GBasisXYZ()
                basis_mp=getMultiParticleBasis(getMultiSiteBasis(basis_sp,s),e)

                # hamiltonian construction
                hamiltonian= MPElectronPerkinsWoelfleHamiltonian(basis_mp, s, U1,U2,J_H)

                # eigensystem
                es=eigensystem(hamiltonian)
                @test length(es[:values])==6
    
                @test (abs.(es[:values]-E1h(U2)*ones(6)).<1e-6*ones(6)) == trues(6)
                
            end
            
            @testset "4e, 1s" begin
                
               # parameters
                e=4
                s=1

                U1=rand(-1.0:10.0)
                J_H=rand(-1.0:10.0)
                U2=U1-2*J_H

                # basis construction
                basis_sp=getT2GBasisXYZ()
                basis_mp=getMultiParticleBasis(getMultiSiteBasis(basis_sp,s),e)

                # hamiltonian construction
                hamiltonian= MPElectronPerkinsWoelfleHamiltonian(basis_mp, s, U1,U2,J_H)

                # eigensystem
                es=eigensystem(hamiltonian) 
                @test length(es[:values])==15
                
                # analytical solutions
                energies=zeros(15)
                energies[1:9]=E2h_1(U2,J_H)*ones(9)
                energies[10:14]=E2h_0(U2,J_H)*ones(5)
                energies[15]=E2h_00(U2,J_H)
                
                # test
                @test (abs.(es[:values]-energies).<1e-6*ones(length(energies)) ) == trues(length(energies))
                
            end
            
            
        end
        
        
        @testset "electron-hole equivalence" begin

            @testset "1s,0h,6e"
                function check_eigensystems(
                        eigensystem1 :: Dict,
                        eigensystem2 :: Dict,
                        index :: Integer
                        ;
                        subtract_GS :: Bool=false,
                        cutoff :: Real = 0.01,
                        digits :: Integer = 3,
                        print_results :: Bool = true
                    )
                    # get the state and the eigenenergy
                    hole_state  = eigensystem1[:vectors][index]
                    hole_energy = eigensystem1[:values][index]

                    electron_state  = eigensystem2[:vectors][index]
                    electron_energy = eigensystem2[:values][index]

                    new_state = 1im.*zeros(length(electron_state))

                    # put minus sign when necessary
                    for i in 1:length(electron_state)
                        occ=basis(eigensystem2[:operator])[i].occupation
                        # add minus sign when sum of occ is odd
                        new_state[i] = (-1)^sum(occ) * electron_state[i]        
                    end

                    if subtract_GS
                        hole_energy -= eigensystem[:values][1]
                        electron_energy -= eigensystem[:values][1]
                    end

                    # compute the overlap; dot already conjugates the first vector
                    ov=dot(hole_state,reverse(new_state))

                    # print the overlap
                    if print_results==true
                        printstyled("The overlap of eigenvector nr. $(index) is:\n", bold=true, color=:light_green)
                        println(round(abs(ov),digits=3))
                    end

                    return (new_state,abs(ov))
                end

                #---------------------------------------------------------------------------------------------------------------
                #---------------------------------------------------------------------------------------------------------------

                # set parameters
                lambda = rand()
                Delta  = rand()
                B_dir=generate_rvos()
                B=rand()

                U=rand()
                J_H=rand()

                P=rand(1:6)


                #---------------------------------------------------------------------------------------------------------------
                #---------------------------------------------------------------------------------------------------------------


                s=1
                p=P
                particle_type=:hole


                ################
                # BASIS
                ################

                # single particle basis for the t2gs
                basis_sp   = getT2GBasisJ()

                # single particle basis for s sites (from single particle basis of t2gs)
                basis_ms = getMultiSiteBasis(basis_sp, s)

                # multi particle basis that occupies single particle 1 site basis with p particles
                basis_mp   = getMultiParticleBasis(basis_ms, p)

                # write interaction basis
                basis_interaction = getMultiParticleBasis(getMultiSiteBasis(getT2GBasisXYZ(), s), p)



                ################
                # HAMILTONIAN
                ################

                # build the sum of all contributing terms (L*S + Delta * Lz^2)
                # Hamiltonian construction

                # compute interaction term in the xyz basis
                if particle_type==:electron
                    interaction = MPElectronPerkinsWoelfleHamiltonian(basis_interaction, 1, 0.0,0.0,0.0)
                else
                    interaction = MPHolePerkinsWoelfleHamiltonian(basis_interaction, 1, 0.0,0.0,0.0)
                end

                # project to chosen basis
                hamiltonian = ProjectorOperator(interaction, basis_mp)

                # SP terms
                K=Pair(:K,0.0)

                hamiltonian += DistortionOperator(basis_mp, 1, 0.0, [0,0,1]; particle_type=particle_type) +
                                SpinOrbitOperator(basis_mp, 1, 0.0; particle_type=particle_type) + 
                                MagneticFieldOperator(basis_mp, 1, 0.0, [0,0,1]; particle_type=particle_type)+
                                K*JzOperator(basis_mp,1; particle_type=particle_type)



                # recalculate the hamiltonian (set matrices to correct values etc.)
                recalculate!(hamiltonian);

                set_parameter!(hamiltonian,      :lambda, lambda, site=:all)
                set_parameter!(hamiltonian,      :Delta,  Delta, site=:all)
                set_parameter!(hamiltonian,      :B,      B, site=:all)
                set_parameter!(hamiltonian,      :B_dir,  B_dir, site=:all)

                set_parameter!(hamiltonian, :U,  U, site=:all)
                set_parameter!(hamiltonian, :J_H,  J_H, site=:all)


                recalculate!(hamiltonian)

                esh=eigensystem(hamiltonian)

                #---------------------------------------------------------------------------------------------------------------
                #---------------------------------------------------------------------------------------------------------------

                s=1
                p=6-P
                particle_type=:electron


                ################
                # BASIS
                ################

                # single particle basis for the t2gs
                basis_sp   = getT2GBasisJ()

                # single particle basis for s sites (from single particle basis of t2gs)
                basis_ms = getMultiSiteBasis(basis_sp, s)

                # multi particle basis that occupies single particle 1 site basis with p particles
                basis_mp   = getMultiParticleBasis(basis_ms, p)

                # write interaction basis
                basis_interaction = getMultiParticleBasis(getMultiSiteBasis(getT2GBasisXYZ(), s), p)



                ################
                # HAMILTONIAN
                ################

                # build the sum of all contributing terms (L*S + Delta * Lz^2)
                # Hamiltonian construction

                # compute interaction term in the xyz basis
                if particle_type==:electron
                    interaction = MPElectronPerkinsWoelfleHamiltonian(basis_interaction, 1, 0.0,0.0,0.0)
                else
                    interaction = MPHolePerkinsWoelfleHamiltonian(basis_interaction, 1, 0.0,0.0,0.0)
                end

                # project to chosen basis
                hamiltonian = ProjectorOperator(interaction, basis_mp)



                # SP terms
                K=Pair(:K,0.0)
            
                hamiltonian += DistortionOperator(basis_mp, 1, 0.0, [0,0,1]; particle_type=particle_type) +
                                SpinOrbitOperator(basis_mp, 1, 0.0; particle_type=particle_type) + 
                                MagneticFieldOperator(basis_mp, 1, 0.0, [0,0,1]; particle_type=particle_type)+
                                K*JzOperator(basis_mp,1; particle_type=particle_type)

                # recalculate the hamiltonian (set matrices to correct values etc.)
                recalculate!(hamiltonian);


                set_parameter!(hamiltonian,      :lambda, lambda, site=:all)
                set_parameter!(hamiltonian,      :Delta,  Delta, site=:all)
                set_parameter!(hamiltonian,      :B,      B, site=:all)
                set_parameter!(hamiltonian,      :B_dir,  B_dir, site=:all)

                set_parameter!(hamiltonian, :U,  U, site=:all)
                set_parameter!(hamiltonian, :J_H,  J_H, site=:all)


                recalculate!(hamiltonian)

                ese=eigensystem(hamiltonian)



                #---------------------------------------------------------------------------------------------------------------
                #---------------------------------------------------------------------------------------------------------------


                # TESTING (test each state; note that the total number of test varies wrt the number of particles, which is randomized)
                sum=0
                for i in 1:length(esh[:vectors])
                    new_state,ov=check_eigensystems(esh,ese,i; print_results=false)

                    sum+=ov

                end
                @test abs(sum-1.0)<1e-6
            
            end # end 1s0h6e
            
        end #end testset e-h equivalence
        
        
    end
    
    
    
end #end physics tests