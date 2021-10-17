# TESTSET FOR ALL ABSTRACT AND CONCRETE TYPES
# mainly test dependencies of abstract types and availability
@testset "Type Definitions" begin
    
    @testset "Basis definition" begin
        
        @testset "Abstract Basis definitions" begin
            @test AbstractBasis <: Any
            @test AbstractBasisState <: Any
            #@test AbstractBasis{BS} <: AbstractArray{BS,1} where {BS<:AbstractBasisState}
        end
        
        @testset "Single-Particle Basis State definitions" begin
            
            @testset "Abstract Single-Paticle Basis State definitions" begin
                @test AbstractSPBasisState <: AbstractBasisState
                @test AbstractSPSSBasisState <: AbstractSPBasisState
            end
            @testset "Concrete Single-Paticle Basis State definitions" begin
                
                @test SPBasis <: AbstractBasis
                
                @testset "Concrete t2g basis state types" begin
                    @test BasisStateA1G <: AbstractSPSSBasisState
                    @test BasisStateJ <: AbstractSPSSBasisState
                    @test BasisStateLS <: AbstractSPSSBasisState
                    @test BasisStateXYZ <: AbstractSPSSBasisState 
                end
                
                @testset "Concrete composite basis state types" begin
                    @test SPSSCompositeBasisState <: AbstractSPSSBasisState
                    @test SPMSCompositeBasisState <: AbstractSPBasisState
                end
                
                @testset "Concrete multi-site basis state types" begin
                    @test SPMSBasisState <: AbstractSPBasisState
                end
                
                
            
            end
        end
            
        
        @testset "Multi-Particle Basis State definitions" begin
        end
        
    end
    
    
    @testset "Abstract Operator definition" begin
            
#             @test AbstractOperator{B <: AbstractBasis{BS} where {BS <: AbstractBasisState}}
            @test AbstractOperator <: Any
            @test AbstractSPOperator <: AbstractOperator
            @test AbstractMPOperator <: AbstractOperator
    end
    
    @testset "Abstract Spectrum and Transition definitions" begin
       
        @test AbstractTransition <: Any
        @test AbstractSpectrum <: Any
        
    end
        
    
# end the testset    
end