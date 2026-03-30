import streamlit as st
import pandas as pd

# ==========================================
# 1. 网页全局设置
# ==========================================
st.set_page_config(page_title="双抗智能组装工厂 V1.0", page_icon="🧬", layout="wide")

st.title("🧬 全构型双抗智能组装与错配防御工厂 (V1.0)")
st.info("💡 架构师模式：输入两个靶点的可变区 (VH/VL)，一键智能拼接生成多链表达序列。自带轻重链错配 (Mispairing) 风险雷达与柔性 Linker 优化引擎。")

# ==========================================
# 2. 核心物料库：标准恒定区与 Linker 字典
# ==========================================
CONSTANT_REGIONS = {
    "CH1": "ASTKGPSVFPLAPSSKSTSGGTAALGCLVKDYFPEPVTVSWNSGALTSGVHTFPAVLQSSGLYSLSSVVTVPSSSLGTQTYICNVNHKPSNTKVDKKV",
    "Hinge": "EPKSCDKTHTCPPCP",
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
# 3. 核心组装引擎
# ==========================================
def assemble_bispecific(ab_a, ab_b, format_type, linker_type):
    chains = []
    alerts = []
    linker = LINKERS[linker_type]
    
    if format_type == "Tandem scFv (BiTE 构型)":
        seq = f"{ab_a['VL']}{linker}{ab_a['VH']}{linker}{ab_b['VL']}{linker}{ab_b['VH']}"
        chains.append({"链名称": "Chain 1 (Full Tandem scFv)", "序列组成": "VL(A)-L-VH(A)-L-VL(B)-L-VH(B)", "完整氨基酸序列": seq})
        alerts.append("✅ 错配风险评估：极低。由于是单链表达，不存在多链错配问题。")
        alerts.append("⚠️ 生产与 CMC 风险：高。BiTE 构型无 Fc 区，半衰期极短，且极易在细胞中发生聚集沉淀。")

    elif format_type == "Asymmetric scFv-IgG (经典不对称 2+1 构型)":
        hc_a = f"{ab_a['VH']}{CONSTANT_REGIONS['CH1']}{CONSTANT_REGIONS['Hinge']}{CONSTANT_REGIONS['CH2_CH3_Knob']}"
        chains.append({"链名称": "Chain 1 (Heavy A - Knob)", "序列组成": "VH(A)-CH1-Hinge-Fc(Knob)", "完整氨基酸序列": hc_a})
        lc_a = f"{ab_a['VL']}{CONSTANT_REGIONS['CL_Kappa']}"
        chains.append({"链名称": "Chain 2 (Light A)", "序列组成": "VL(A)-CL", "完整氨基酸序列": lc_a})
        scfv_b_fc = f"{ab_b['VL']}{linker}{ab_b['VH']}{CONSTANT_REGIONS['Hinge']}{CONSTANT_REGIONS['CH2_CH3_Hole']}"
        chains.append({"链名称": "Chain 3 (scFv B - Hole)", "序列组成": "VL(B)-L-VH(B)-Hinge-Fc(Hole)", "完整氨基酸序列": scfv_b_fc})
        
        alerts.append("✅ 错配风险评估：低。有效规避了轻重链交叉错配 (Light-chain mispairing)。")
        alerts.append("⚠️ 结构空间风险：中。需注意 scFv(B) 融合在铰链区可能会与对侧的 Fab(A) 产生空间位阻。")

    elif format_type == "Dual-Variable Domain Ig (DVD-Ig 构型)":
        hc_dvd = f"{ab_b['VH']}{linker}{ab_a['VH']}{CONSTANT_REGIONS['CH1']}{CONSTANT_REGIONS['Hinge']}{CONSTANT_REGIONS['CH2_CH3_Knob']}"
        chains.append({"链名称": "Chain 1 (DVD Heavy)", "序列组成": "VH(B)-L-VH(A)-CH1-Hinge-Fc", "完整氨基酸序列": hc_dvd})
        lc_dvd = f"{ab_b['VL']}{linker}{ab_a['VL']}{CONSTANT_REGIONS['CL_Kappa']}"
        chains.append({"链名称": "Chain 2 (DVD Light)", "序列组成": "VL(B)-L-VL(A)-CL", "完整氨基酸序列": lc_dvd})
        
        alerts.append("✅ 错配风险评估：极低。对称的四价结构，只需表达两条链。")
        alerts.append("🚨 靶点结合风险：极高。内侧的抗原结合域 (Paratope A) 可能被严重遮挡，导致亲和力断崖式下降。")

    elif format_type == "Classic IgG-like Bispecific (无保护经典组装)":
        hc_a = f"{ab_a['VH']}{CONSTANT_REGIONS['CH1']}{CONSTANT_REGIONS['Hinge']}{CONSTANT_REGIONS['CH2_CH3_Knob']}"
        lc_a = f"{ab_a['VL']}{CONSTANT_REGIONS['CL_Kappa']}"
        hc_b = f"{ab_b['VH']}{CONSTANT_REGIONS['CH1']}{CONSTANT_REGIONS['Hinge']}{CONSTANT_REGIONS['CH2_CH3_Hole']}"
        lc_b = f"{ab_b['VL']}{CONSTANT_REGIONS['CL_Kappa']}"
        chains = [
            {"链名称": "Chain 1 (Heavy A - Knob)", "序列组成": "VH(A)-CH1-Hinge-Fc(Knob)", "完整氨基酸序列": hc_a},
            {"链名称": "Chain 2 (Light A)", "序列组成": "VL(A)-CL", "完整氨基酸序列": lc_a},
            {"链名称": "Chain 3 (Heavy B - Hole)", "序列组成": "VH(B)-CH1-Hinge-Fc(Hole)", "完整氨基酸序列": hc_b},
            {"链名称": "Chain 4 (Light B)", "序列组成": "VL(B)-CL", "完整氨基酸序列": lc_b}
        ]
        alerts.append("🚨 错配风险评估：致命级！两根不同的轻链会在细胞内随机与两根重链结合，产物中会有 50% 以上的杂质！")
        alerts.append("💡 架构师建议：必须引入 CrossMab 技术，或者使用 Common Light Chain（共用轻链）策略进行底层重构。")

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
st.markdown("### 🏗️ 第二步：选择空间拓扑与组装策略")

col3, col4 = st.columns(2)
with col3:
    format_choice = st.selectbox("选择双抗结构拓扑 (Format):", [
        "Asymmetric scFv-IgG (经典不对称 2+1 构型)", 
        "Tandem scFv (BiTE 构型)", 
        "Dual-Variable Domain Ig (DVD-Ig 构型)",
        "Classic IgG-like Bispecific (无保护经典组装)"
    ])
with col4:
    linker_choice = st.selectbox("选择内部柔性接头 (Linker):", list(LINKERS.keys()))

if st.button("🚀 启动智能拼装工厂", type="primary", use_container_width=True):
    ab_a = {"VH": vh_a.replace(".", "").replace(" ", "").upper(), "VL": vl_a.replace(".", "").replace(" ", "").upper()}
    ab_b = {"VH": vh_b.replace(".", "").replace(" ", "").upper(), "VL": vl_b.replace(".", "").replace(" ", "").upper()}
    
    if len(ab_a["VH"]) < 50 or len(ab_b["VH"]) < 50:
        st.error("❌ 序列异常：请确保输入的 VH 和 VL 长度有效 (至少 50 aa 以上)。")
    else:
        chains, alerts = assemble_bispecific(ab_a, ab_b, format_choice, linker_choice)
        
        st.markdown("---")
        st.markdown(f"### 📦 组装产物清单: {format_choice}")
        
        for alert in alerts:
            if "✅" in alert:
                st.success(alert)
            elif "⚠️" in alert or "💡" in alert:
                st.warning(alert)
            elif "🚨" in alert:
                st.error(alert)
                
        df_chains = pd.DataFrame(chains)
        st.dataframe(df_chains, use_container_width=True)
        
        fasta_output = ""
        for chain in chains:
            fasta_output += f">{chain['链名称']} | {chain['序列组成']}\n{chain['完整氨基酸序列']}\n\n"
            
        st.markdown("#### 📋 质粒合成级 FASTA 序列")
        st.text_area("可直接全选复制，发送给基因合成供应商：", value=fasta_output, height=250)
