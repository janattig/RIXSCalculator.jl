# TESTSET FOR ALL ABSTRACT AND CONCRETE TYPES
# mainly test dependencies of abstract types and availability
@testset "Type Definitions" begin
    
    @testset "Basis definition" begin
        
        @testset "Abstract Basis definitions" begin
            @testset "Abstract Single-Particle Basis definitions" begin
                @test AbstractBasis <: Any
                @test AbstractBasisState <: Any

                @test AbstractSPBasisState <: AbstractBasisState
                @test AbstractSPSSBasisState <: AbstractSPBasisState

                @test SPBasis <: AbstractBasis
            end
        end
        
        @testset "Concrete Basis State definitions" begin
            
            @testset "Single-Particle Basis State definitions" begin
                
                @testset "t2g basis state types" begin
                    @test BasisStateA1G <: AbstractSPSSBasisState
                    @test BasisStateJ <: AbstractSPSSBasisState
                    @test BasisStateLS <: AbstractSPSSBasisState
                    @test BasisStateXYZ <: AbstractSPSSBasisState 
                end
                
                @testset "composite basis state types" begin
                    @test SPSSCompositeBasisState <: AbstractSPSSBasisState
                    @test SPMSCompositeBasisState <: AbstractSPBasisState
                end
                
                @testset "multi-site basis state types" begin
                    @test SPMSBasisState <: AbstractSPBasisState
                end
            end
            
            
            @testset "Multi-Particle Basis State definitions" begin
                
                @test MPBasisState <: AbstractBasisState
                @test MPBasis <: AbstractBasis
                
            end
            
        end #end concrete basis state definitions
            
        
    end
  
#######################################################################################
    
    
    @testset "Operator definition" begin
        
        @testset "Abstract Operator definitions" begin
            @test AbstractOperator <: Any
            @test AbstractSPOperator <: AbstractOperator
            @test AbstractMPOperator <: AbstractOperator
        end

            
    end
    
    
#######################################################################################   
    
    @testset "Abstract Spectrum and Transition definitions" begin
       
        @test AbstractTransition <: Any
        @test AbstractSpectrum <: Any
        
    end
        
    
# end the testset    
end