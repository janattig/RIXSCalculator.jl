# TESTSET FOR ALL ABSTRACT TYPES
# mainly test dependencies of abstract types and availability
@testset "Abstract Type Definitions" begin
    
    @testset "Abstract Basis definition" begin
        
        @test AbstractBasis{BS} <: AbstractArray{BS,1} where {BS}
        @test_nowarn AbstractBasisState <: Any
        
        @testset "Single-Particle Basis State definitions" begin
            
            @test AbstractSPBasisState <: AbstractBasisState
            @test AbstractSPSSBasisState <: AbstractSPBasisState
            
        end
        
        @testset "Multi-Particle Basis State definitions" begin
        end
        
    end
    
    
    @testset "Abstract Operator definition" begin
            
#             @test AbstractOperator{B <: AbstractBasis{BS} where {BS <: AbstractBasisState}}
        
            @test AbstractSPOperator <: AbstractOperator
            @test AbstractMPOperator <: AbstractOperator
    end
    
    @testset "Abstract Spectrum and Transition definitions" begin
       
        @test AbstractTransition
        @test AbstractSpectrum
        
    end
        
    
# end the testset    
end