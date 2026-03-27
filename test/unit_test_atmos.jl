@testset "TASOPT.atmospheric model" begin
    @testset "Sea level reference" begin
        T, p, ρ, a, μ = TASOPT.atmos(0.0)

        @test T ≈ 288.15 rtol=1e-5
        @test p ≈ 1.0132e5 rtol=1e-4
        @test ρ ≈ 1.2250 rtol=1e-5
        @test a ≈ 340.294 rtol=1e-5
        @test μ ≈ 1.7894e-5 rtol=1e-5
    end

    @testset "Layer interfaces" begin
        for h_int in [11.0, 20.0, 32.0]
            Tm, pm, ρm, am, μm = TASOPT.atmos(h_int - 1e-6)
            Tp, pp, ρp, ap, μp = TASOPT.atmos(h_int + 1e-6)

            @test Tm ≈ Tp atol=1e-4
            @test pm ≈ pp rtol=1e-4
            @test ρm ≈ ρp rtol=1e-4
        end
    end

    @testset "Lower stratosphere reference" begin
        T, p, ρ, a, μ = TASOPT.atmos(20.0)

        @test T ≈ 216.65 rtol=1e-6
        @test p ≈ 5474.9 rtol=1e-3
        @test ρ ≈ 0.0880 rtol=1e-3
        @test a ≈ 295.1 rtol=5e-4
    end

    @testset "Stratospheric warming" begin
        T20, p20, ρ20, a20, μ20 = TASOPT.atmos(20.0)
        T32, p32, ρ32, a32, μ32 = TASOPT.atmos(32.0)
        T47, p47, ρ47, a47, μ47 = TASOPT.atmos(47.0)

        @test T32 > T20
        @test T47 > T32
        @test p32 < p20
        @test p47 < p32
        @test ρ32 < ρ20
        @test ρ47 < ρ32
    end

    @testset "Nonstandard temperature" begin
        Tstd, pstd, ρstd, astd, μstd = TASOPT.atmos(20.0, 0.0)
        Thot, phot, ρhot, ahot, μhot = TASOPT.atmos(20.0, 10.0)

        @test Thot ≈ Tstd + 10.0
        @test phot ≈ pstd rtol=1e-12
        @test ρhot < ρstd
        @test ahot > astd
        @test μhot > μstd
    end

    @testset "Monotonic pressure decrease" begin
        hs = [0.0, 5.0, 11.0, 15.0, 20.0, 25.0, 32.0, 40.0, 47.0]
        ps = [TASOPT.atmos(h)[2] for h in hs]

        @test all(ps[i] > ps[i+1] for i in 1:length(ps)-1)
    end

    @testset "Altitude from density" begin
        ρ0 = 1.2250
        ρ10 = 0.41271
        ρ20 = 0.08803
        ρ30 = 0.01801
        ρ40 = 0.003851
        
        @test TASOPT.find_altitude_from_density(ρ0) ≈ 0.0 atol=1e-3
        @test TASOPT.find_altitude_from_density(ρ10) ≈ 10.0 atol=1e-3
        @test TASOPT.find_altitude_from_density(ρ20) ≈ 20.0 atol=1e-3
        @test TASOPT.find_altitude_from_density(ρ30) ≈ 30.0 atol=1e-3
        @test TASOPT.find_altitude_from_density(ρ40) ≈ 40.0 atol=1e-3
    end
end