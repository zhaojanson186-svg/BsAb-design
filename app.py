import streamlit as st
import pandas as pd
from Bio.SeqUtils.ProtParam import ProteinAnalysis

# ==========================================
# 1. 网页全局设置
# ==========================================
st.set_page_config(page_title="双抗智能组装工厂 V3.0", page_icon="🧬", layout="wide")

st.title("🧬 全构型双抗智能组装与 CMC 预测工厂 (V3.0)")
st.info("💡 终极智造模式：内置 7 大拓扑构型，新增【单链 pI/MW 物理化学属性即时测算】与【一键 LALA-PG 毒性沉默】工程插件。")

# ==========================================
# 2. 核心物料库：标准恒定区与 Linker 字典
# ==========================================
CONSTANT_REGIONS = {
    "CH1": "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV",
    "Hinge": "EPKSCDKTHTCPPCP",
    "CH2_CH3": "APELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLTCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK", 
    "CH2_CH3_Knob": "APELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLWCLVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLYSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK", 
    "CH2_CH3_Hole": "APELLGGPSVFLFPPKPKDTLMISRTPEVTCVVVDVSHEDPEVKFNWYVDGVEVHNAKTKPREEQYNSTYRVVSVLTVLHQDWLNGKEYKCKVSNKALPAPIEKTISKAKGQPREPQVYTLPPSREEMTKNQVSLSCAVKGFYPSDIAVEWESNGQPENNYKTTPPVLDSDGSFFLVSKLTVDKSRWQQGNVFSCSVMHEALHNHYTQKSLSLSPGK", 
    "CL_Kappa": "RTVAAPSVFIFPPSDEQLKSGTASVVCLLNNFYPREAKVQWKVDNALQSGNSQESVTEQDSKDSTYSLSSTLTLSKADYEKHKVYACEVTHQGLSSPVTKSFNRGEC"
}

LINKERS = {
    "GS15 (经典的 3x G4S)": "GGGGSGGGGSGGGGS",
    "GS20 (柔性更强的 4x G4S)": "GGGGSGGGGSGGGGSGGGGS",
    "GS5 (极短刚性 1x G4S)": "GGGGS",
    "ASTK (天然构象倾向)": "ASTKGPSVFPLAP"
}

# ==========================================
# 3. 核心计算与组装引擎
# ==========================================
def calculate_cmc_props(seq):
    """利用 Biopython 计算分子量和等电点"""
    try:
        cleaned_seq = seq.replace("X", "").replace("-", "")
        pa = ProteinAnalysis(cleaned_seq)
        mw = round(pa.molecular_weight() / 1000, 2) # kDa
        pi = round(pa.isoelectric_point(), 2)
        return mw, pi
    except:
        return 0.0, 0.0

def apply_lala_pg(seq):
    """精准识别并替换 L234A, L235A 和 P329G"""
    modified = seq.replace("PELLGG", "PEAAGG") # LALA
    modified = modified.replace("KALPAPI", "KALGAPI") # P329G
    return modified

def assemble_bispecific(ab_a, ab_b, format_type, linker_type, use_lalapg):
    chains = []
    alerts = []
    linker = LINKERS[linker_type]
    
    # 动态获取常域，如果开启了静默，则动态替换序列
    fc_wt = apply_lala_pg(CONSTANT_REGIONS["CH2_CH3"]) if use_lalapg else CONSTANT_REGIONS["CH2_CH3"]
    fc_knob = apply_lala_pg(CONSTANT_REGIONS["CH2_CH3_Knob"]) if use_lalapg else CONSTANT_REGIONS["CH2_CH3_Knob"]
    fc_hole = apply_lala_pg(CONSTANT_REGIONS["CH2_CH3_Hole"]) if use_lalapg else CONSTANT_REGIONS["CH2_CH3_Hole"]
    
    # 智能识别 Common Light Chain (cLC)
    is_clc = ab_a['VL'] == ab_b['VL']
    if is_clc and format_type != "Tandem scFv (BiTE 构型)":
        alerts.append("💡 智能识别：检测到输入的抗体 A 和 B 使用了【完全相同的轻链 (cLC)】！这从根本上消除了轻链错配问题，极大地简化了下游纯化与放大工艺。")

    if format_type == "Tandem scFv (BiTE 构型)":
        seq = f"{ab_a['VL']}{linker}{ab_a['VH']}{linker}{ab_b['VL']}{linker}{ab_b['VH']}"
        chains.append({"链名称": "Chain 1 (Full Tandem scFv)", "序列组成": "VL(A)-L-VH(A)-L-VL(B)-L-VH(B)", "完整氨基酸序列": seq})
        alerts.append("✅ 错配风险评估：极低。单链表达，不存在多链错配。")
        alerts.append("⚠️ CMC 风险：高。无 Fc 区导致半衰期极短，且极易在细胞表达中发生聚集。")

    elif format_type == "Asymmetric scFv-IgG (经典不对称 2+1 构型)":
        hc_a = f"{ab_a['VH']}{CONSTANT_REGIONS['CH1']}{CONSTANT_REGIONS['Hinge']}{fc_knob}"
        chains.append({"链名称": "Chain 1 (Heavy A - Knob)", "序列组成": "VH(A)-CH1-Hinge-Fc(Knob)", "完整氨基酸序列": hc_a})
        lc_a = f"{ab_a['VL']}{CONSTANT_REGIONS['CL_Kappa']}"
        chains.append({"链名称": "Chain 2 (Light A)", "序列组成": "VL(A)-CL", "完整氨基酸序列": lc_a})
        scfv_b_fc = f"{ab_b['VL']}{linker}{ab_b['VH']}{CONSTANT_REGIONS['Hinge']}{fc_hole}"
        chains.append({"链名称": "Chain 3 (scFv B - Hole)", "序列组成": "VL(B)-L-VH(B)-Hinge-Fc(Hole)", "完整氨基酸序列": scfv_b_fc})
        
        alerts.append("✅ 错配风险评估：低。一端为天然 Fab，一端为单链，有效规避了交叉错配。")

    elif format_type == "CrossMab (CH1/CL Swap 经典 1+1 构型)":
        hc_a = f"{ab_a['VH']}{CONSTANT_REGIONS['CH1']}{CONSTANT_REGIONS['Hinge']}{fc_knob}"
        chains.append({"链名称": "Chain 1 (Heavy A - Knob)", "序列组成": "VH(A)-CH1-Hinge-Fc(Knob)", "完整氨基酸序列": hc_a})
        lc_a = f"{ab_a['VL']}{CONSTANT_REGIONS['CL_Kappa']}"
        chains.append({"链名称": "Chain 2 (Light A)", "序列组成": "VL(A)-CL", "完整氨基酸序列": lc_a})
        
        hc_b_swap = f"{ab_b['VH']}{CONSTANT_REGIONS['CL_Kappa']}{CONSTANT_REGIONS['Hinge']}{fc_hole}"
        chains.append({"链名称": "Chain 3 (Heavy B Swap - Hole)", "序列组成": "VH(B)-CL-Hinge-Fc(Hole)", "完整氨基酸序列": hc_b_swap})
        lc_b_swap = f"{ab_b['VL']}{CONSTANT_REGIONS['CH1']}"
        chains.append({"链名称": "Chain 4 (Light B Swap)", "序列组成": "VL(B)-CH1", "完整氨基酸序列": lc_b_swap})
        
        alerts.append("✅ 错配风险评估：极低。利用界面的拓扑排斥，完美解决了轻链错配问题 (Roche 核心专利)。")

    elif format_type == "IgG-scFv (C端融合 四价 2+2 构型)":
        hc_igg_scfv = f"{ab_a['VH']}{CONSTANT_REGIONS['CH1']}{CONSTANT_REGIONS['Hinge']}{fc_wt}{linker}{ab_b['VL']}{linker}{ab_b['VH']}"
        chains.append({"链名称": "Chain 1 (Heavy A + scFv B)", "序列组成": "VH(A)-CH1-Hinge-Fc-L-VL(B)-L-VH(B)", "完整氨基酸序列": hc_igg_scfv})
        lc_a = f"{ab_a['VL']}{CONSTANT_REGIONS['CL_Kappa']}"
        chains.append({"链名称": "Chain 2 (Light A)", "序列组成": "VL(A)-CL", "完整氨基酸序列": lc_a})
        
        alerts.append("✅ 错配风险评估：极低。对称的同源二聚体，细胞只需表达 2 条链。")
        alerts.append("🚨 降解风险：极高！C端悬挂的 scFv 非常容易在细胞表达或体内循环时发生断裂 (Clipping)。")

    elif format_type == "Asymmetric scFv-Fc (小分子量 1+1 构型)":
        scfv_a_fc = f"{ab_a['VL']}{linker}{ab_a['VH']}{CONSTANT_REGIONS['Hinge']}{fc_knob}"
        chains.append({"链名称": "Chain 1 (scFv A - Knob)", "序列组成": "VL(A)-L-VH(A)-Hinge-Fc(Knob)", "完整氨基酸序列": scfv_a_fc})
        scfv_b_fc = f"{ab_b['VL']}{linker}{ab_b['VH']}{CONSTANT_REGIONS['Hinge']}{fc_hole}"
        chains.append({"链名称": "Chain 2 (scFv B - Hole)", "序列组成": "VL(B)-L-VH(B)-Hinge-Fc(Hole)", "完整氨基酸序列": scfv_b_fc})
        
        alerts.append("✅ 错配风险评估：极低。全单链抗体组装，完全无轻链干扰。")

    elif format_type == "Dual-Variable Domain Ig (DVD-Ig 构型)":
        hc_dvd = f"{ab_b['VH']}{linker}{ab_a['VH']}{CONSTANT_REGIONS['CH1']}{CONSTANT_REGIONS['Hinge']}{fc_wt}" 
        chains.append({"链名称": "Chain 1 (DVD Heavy)", "序列组成": "VH(B)-L-VH(A)-CH1-Hinge-Fc", "完整氨基酸序列": hc_dvd})
        lc_dvd = f"{ab_b['VL']}{linker}{ab_a['VL']}{CONSTANT_REGIONS['CL_Kappa']}"
        chains.append({"链名称": "Chain 2 (DVD Light)", "序列组成": "VL(B)-L-VL(A)-CL", "完整氨基酸序列": lc_dvd})
        
        alerts.append("✅ 错配风险评估：极低。对称四价结构。")
        alerts.append("🚨 靶点结合风险：极高。内侧抗原结合域极易产生严重的空间遮挡 (Steric Hindrance)。")

    elif format_type == "Classic IgG-like Bispecific (无保护经典组装)":
        hc_a = f"{ab_a['VH']}{CONSTANT_REGIONS['CH1']}{CONSTANT_REGIONS['Hinge']}{fc_knob}"
        lc_a = f"{ab_a['VL']}{CONSTANT_REGIONS['CL_Kappa']}"
        hc_b = f"{ab_b['VH']}{CONSTANT_REGIONS['CH1']}{CONSTANT_REGIONS['Hinge']}{fc_hole}"
        lc_b = f"{ab_b['VL']}{CONSTANT_REGIONS['CL_Kappa']}"
        chains = [
            {"链名称": "Chain 1 (Heavy A - Knob)", "序列组成": "VH(A)-CH1-Hinge-Fc(Knob)", "完整氨基酸序列": hc_a},
            {"链名称": "Chain 2 (Light A)", "序列组成": "VL(A)-CL", "完整氨基酸序列": lc_a},
            {"链名称": "Chain 3 (Heavy B - Hole)", "序列组成": "VH(B)-CH1-Hinge-Fc(Hole)", "完整氨基酸序列": hc_b},
            {"链名称": "Chain 4 (Light B)", "序列组成": "VL(B)-CL", "完整氨基酸序列": lc_b}
        ]
        if not is_clc:
            alerts.append("🚨 错配风险评估：致命级！两根不同的轻链会在细胞内随机与两根重链结合，产物中会有 50% 以上的杂质！建议切换至 CrossMab。")

    # 追加计算每一条链的 CMC 属性
    for chain in chains:
        mw, pi = calculate_cmc_props(chain["完整氨基酸序列"])
        chain["MW 分子量 (kDa)"] = mw
        chain["pI 等电点"] = pi

    return chains, alerts

# ==========================================
# 4. 交互界面设计
# ==========================================
st.markdown("### 🧰 第一步：录入抗体 Fv 构件")
st.caption("请分别输入靶向端和效应端的可变区序列。工具将自动为您挂载标准 IgG1 骨架。")

col1, col2 = st.columns(2)
with col1:
    st.markdown("#### 🎯 抗体 A (如 Tumor Target)")
    vh_a = st.text_area("VH (Heavy Chain Variable):", value="EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAK", height=100, key="vha")
    vl_a = st.text_area("VL (Light Chain Variable):", value="DIQMTQSPSSLSASVGDRVTITCRASQGISNYLAWYQQKPGKAPKLLIYAASTLQSGVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQLNSYPLT", height=100, key="vla")

with col2:
    st.markdown("#### ⚔️ 抗体 B (如 CD3 Effector)")
    vh_b = st.text_area("VH (Heavy Chain Variable):", value="EVQLVESGGGLVQPGGSLKLSCAASGFTFNKYAMNWVRQAPGKGLEWVARIRSKYNNYATYYADSVKDRFTISRDDSKNTAYLQMNNLKTEDTAVYYCVR", height=100, key="vhb")
    vl_b = st.text_area("VL (Light Chain Variable):", value="QTVVTQEPSLTVSPGGTVTLTCGSSTGAVTSGNYPNWVQQKPGQAPRGLIGGTKFLAPGTPARFSGSLLGGKAALTLSGVQPEDEAEYYCVLWYSNRW", height=100, key="vlb")

st.markdown("---")
st.markdown("### 🏗️ 第二步：选择空间拓扑与高级工程插件")

col3, col4, col5 = st.columns([2, 1, 1])
with col3:
    format_choice = st.selectbox("选择双抗结构拓扑 (Format):", [
        "CrossMab (CH1/CL Swap 经典 1+1 构型)",
        "Asymmetric scFv-IgG (经典不对称 2+1 构型)", 
        "IgG-scFv (C端融合 四价 2+2 构型)",
        "Asymmetric scFv-Fc (小分子量 1+1 构型)",
        "Tandem scFv (BiTE 构型)", 
        "Dual-Variable Domain Ig (DVD-Ig 构型)",
        "Classic IgG-like Bispecific (无保护经典组装)"
    ])
with col4:
    linker_choice = st.selectbox("内部柔性接头 (Linker):", list(LINKERS.keys()))
with col5:
    st.markdown("<div style='margin-top: 30px;'></div>", unsafe_allow_html=True)
    use_lalapg = st.checkbox("🛡️ 引入 LALA-PG 去毒性突变", value=True, help="勾选后，底层的 IgG1 Fc 区将自动替换为带有 L234A, L235A, P329G 的静默版本，大幅降低 ADCC 毒性。")

if st.button("🚀 启动智能拼装工厂", type="primary", use_container_width=True):
    ab_a = {"VH": vh_a.replace(".", "").replace(" ", "").upper(), "VL": vl_a.replace(".", "").replace(" ", "").upper()}
    ab_b = {"VH": vh_b.replace(".", "").replace(" ", "").upper(), "VL": vl_b.replace(".", "").replace(" ", "").upper()}
    
    if len(ab_a["VH"]) < 50 or len(ab_b["VH"]) < 50:
        st.error("❌ 序列异常：请确保输入的 VH 和 VL 长度有效 (至少 50 aa 以上)。")
    else:
        chains, alerts = assemble_bispecific(ab_a, ab_b, format_choice, linker_choice, use_lalapg)
        
        st.markdown("---")
        st.markdown(f"### 📦 组装产物与 CMC 评估表: `{format_choice}`")
        
        for alert in alerts:
            if "✅" in alert:
                st.success(alert)
            elif "⚠️" in alert or "💡" in alert:
                st.warning(alert)
            elif "🚨" in alert:
                st.error(alert)
                
        # 整理 DataFrame 列顺序，优先展示关键 CMC 数据
        df_chains = pd.DataFrame(chains)
        cols = ['链名称', 'MW 分子量 (kDa)', 'pI 等电点', '序列组成', '完整氨基酸序列']
        df_chains = df_chains[cols]
        
        # 凸显 pI 列以供下游纯化参考
        st.dataframe(df_chains.style.format({'MW 分子量 (kDa)': '{:.2f}', 'pI 等电点': '{:.2f}'}).background_gradient(subset=['pI 等电点'], cmap='viridis'), use_container_width=True)
        
        fasta_output = ""
        for chain in chains:
            fasta_output += f">{chain['链名称']} | {chain['序列组成']} | MW:{chain['MW 分子量 (kDa)']}kDa | pI:{chain['pI 等电点']}\n{chain['完整氨基酸序列']}\n\n"
            
        st.markdown("#### 📋 质粒合成级 FASTA 序列 (含 CMC 标签)")
        st.text_area("可直接全选复制，发送给基因合成供应商：", value=fasta_output, height=250)
