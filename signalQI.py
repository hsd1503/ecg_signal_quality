import numpy as np

def preScreening(sig, rpos):
    # 信号质量粗筛

    # 判断R波数量,若小于6个则判断信号质量差,可根据信号长度调整该参数
    if len(rpos) < 6:
        return False
    
    # 计算最大最小幅值间距,进而初步判断信号质量,若大于5mV则判断信号质量差
    ampl = np.abs(np.max(sig) - np.min(sig))
    if ampl > 5:
       return False
    
    # 核查R波和S波幅值,若幅值过小,判断信号质量差
    tmp_sig = np.abs(sig)
    tmp_sig = tmp_sig[tmp_sig > 0.1]
    if len(tmp_sig) < 5:
        return False
    
    # 比较最大RR间期和平均RR间期,进而初步判断信号质量
    rr_intervals = np.diff(rpos)
    maxRR = np.max(rr_intervals)
    meanRR = np.mean(rr_intervals)
    if maxRR > meanRR * 3:
        return False
    
    return rpos

def cal_corr_coeff_lst(sig, rpos, fs=121):
    # 计算相关系数
    beat_seg = []

    for idx, r_p in enumerate(rpos):
        # 截取心搏
        # 判断第一个R波位置之前是否能够取值,若不能则丢弃
        if idx == 0:
            if r_p > int(0.4*fs):
                tmp_seg = sig[r_p - int(0.4*fs):r_p + int(0.4*fs)]
                beat_seg.append(tmp_seg)
        # 判断最后一个R波位置之后是否能够取值,若不能则丢弃
        elif idx == len(rpos) - 1:
            if (len(sig) - r_p) > int(0.4*fs):
                tmp_seg = sig[r_p - int(0.4*fs):r_p + int(0.4*fs)]
                beat_seg.append(tmp_seg)
        # 截取心搏
        else:
            tmp_seg = sig[r_p - int(0.4*fs):r_p + int(0.4*fs)]
            beat_seg.append(tmp_seg)
    
    # 计算平均心搏
    beat_seg = np.array(beat_seg)
    template_qrs = np.mean(beat_seg, axis=0)
    template_qrs = template_qrs - np.mean(template_qrs)

    # 计算心搏相关系数
    coeff_lst = []
    for seg in beat_seg:
        seg = seg - np.mean(seg)
        coeff = np.corrcoef(seg, template_qrs)[0, 1]
        coeff_lst.append(coeff)
    return coeff_lst

def check_coeff(coeffs):
    # 核查相关系数
    coeffs = np.array(coeffs)
    count = 0

    # 遍历相关系数,通过计数判断相关系数高低进而推算信号质量
    while len(coeffs) > 0:
        coeffs = np.abs(coeffs - coeffs[0])
        coeffs = coeffs[coeffs > 0.1]
        count += 1
    return count

def sqi(sig, rpos, fs=121):
    # 粗筛信号质量
    preRes = preScreening(sig, rpos)

    if preRes:
        rpos = rpos
    else:
        return False, float(0)
    
    # 计算心搏相关系数
    coeff_lst = cal_corr_coeff_lst(sig, rpos, fs)
    template_nums = check_coeff(coeff_lst)
    coeff = float(np.mean(coeff_lst))

    # 判断信号质量
    if coeff > 0.8 and template_nums < 3:
        return True, coeff
    else:
        return False, coeff

with open('ecg.txt', 'r') as fr:
    data = fr.read()
data = [float(x) for x in data.split(",")] 

rposlist = [98, 186, 273, 358, 445, 527, 610, 693, 778, 875, 961, 1049, 1134, 1217, 1301, 1385, 1470, 1554, 1636, 1719, 1802, 1901, 1983, 2069, 2150, 2232]

ecgsqi = sqi(data, rposlist)
print(ecgsqi)