# Rplot

# 1. 概述

本项目应用计算生物学算法进行蛋白质翻译后修饰质谱数据位点匹配以进行新型质谱图绘制并封装成R包。

# 2. 输入
input_table structure like this：
Raw file	Scan number	Sequence	Modifications	mod	Fixed_mod	Gene Names	label
test/20101227_Velos1_TaGe_SA_GAMG_101230100451.mzML	10774	GGGGNFGPGPGSNFR			C:Carbamidomethyl		1
test/20101227_Velos1_TaGe_SA_GAMG_101230100451.mzML	10874	GGGGNFGPGPGSNFR	Oxid_11	F6:ox	C:Carbamidomethyl		1
test/20101227_Velos1_TaGe_SA_GAMG_101230100451.mzML	10874	GGGGNFGPGPGSNFR	Oxid_11	P8:ox	C:Carbamidomethyl		1
test/20101227_Velos1_TaGe_SA_GAMG_101230100451.mzML	10874	GGGGNFGPGPGSNFR	Oxid_11	P10:ox	C:Carbamidomethyl		1

# 3. 输出
