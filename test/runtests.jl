# ~/~ begin <<docs/index.md#test/runtests.jl>>[init]
using Test
using MindTheGap.Stencil
using MindTheGap.Utility

@testset "Mind The Gap tests" begin
    # ~/~ begin <<docs/stencils.md#spec>>[init]
    @testset "offset_value" begin
        @test CartesianIndex(1, 1) == offset_index(Reflected{2}, (3, 3), CartesianIndex(1, 1), CartesianIndex(0, 0))
    end
    # ~/~ end
end
# ~/~ end