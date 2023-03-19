using LinearAlgebra
#########################
# jd = Dates.datetime2rata(DateTime(2010,1,1,0))
#############
# Explan. Suppl.中的系数
_sc = [270.434164, 13.1763965268, -0.0000850, 0.000000039]
_hc = [279.696678, 0.9856473354, 0.00002267, 0.000000000]
_pc = [334.329556, 0.1114040803, -0.0007739, -0.00000026]
_npc = [-259.183275, 0.0529539222, -0.0001557, -0.000000050]
# Explan. Suppl.中的第一个系数是281.220833，但Foreman中是44.
_ppc = [281.220844, 0.0000470684, 0.0000339, 0.000000070]
_coefs = vcat(_sc', _hc', _pc', _npc', _ppc')

function ut_astron(jd)
    jd = vec(jd)  # 确保jd是一维数组
    # 将时代移至1899年12月31日中午：
    # daten = 693961.500000000，Matlab datenum版本
    daten = 693595.5  # Python时代比Matlab晚366天
    d = jd .- daten
    D = d ./ 10000
    args = vcat(ones(length(jd)), d, D .* D, D .^ 3)
    astro = mod.((_coefs * args) ./ 360, 1)  # 月时：太阳日的小数部分加上到太阳经度的时角，减去月亮经度
    tau = mod.(jd, 1) .+ astro[2, :] .- astro[1, :]
    astro = vcat(tau, astro)  # 将tau和astro垂直堆叠
    # 求导数（多项式）
    dargs = vcat(zeros(length(jd)), ones(length(jd)), 2.0e-4 .* D, 3.0e-4 .* D .* D)
    ader = (_coefs * dargs) ./ 360
    dtau = 1.0 .+ ader[2, :] .- ader[1, :]
    ader = vcat(dtau, ader)  # 将dtau和ader垂直堆叠
    return astro, ader
end