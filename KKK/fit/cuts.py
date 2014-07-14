# Custom cuts
cuts_ = " DTFchi2ndof > 0"
cuts_ += "&& DTFchi2ndof < 5"
cuts_ += "&& DTFctau > 0.25"
cuts_ += "&& vchi2_b < 20"
cuts_ += "&& pt_kaon[0] > 0.6 && pt_kaon[1] > 0.6 && pt_kaon[2] > 0.6"
cuts_ += "&& m_jpsi    > 3.020 && m_jpsi    < 3.135"
cuts_ += "&& minann_K  > 0.3"
cuts_ += "&& MIPCHI2DV_k1 > 12. && MIPCHI2DV_k2 > 12. && MIPCHI2DV_k3 > 12."
cuts_ += "&& ann_kaon_PI[0] < 0.3 && ann_kaon_PI[1] < 0.3 &&ann_kaon_PI[2] < 0.3"
# Trigger
cuts_ += "&& ((psi_l0tos & 2) == 2)"
cuts_ += "&& ((psi_l1tos & 2) == 2)"
cuts_ += "&& ((psi_l2tos & 2) == 2)"



mctrue = " && mcTrueB && mcTruePsi && mcTrueK1 && mcTrueK2 && mcTrueK3 && mcTrueMu1 && mcTrueMu2"

cuts_mc = cuts_ + mctrue

# old Cuts
# cuts_ = " DTFchi2ndof > 0"
# cuts_ += "&& DTFchi2ndof < 5"
# cuts_ += "&& DTFctau > 0.2"
# cuts_ += "&& pt_kaon[0] > 0.6 && pt_kaon[1] > 0.6 && pt_kaon[2] > 0.6 "
# cuts_ += "&& m_jpsi    > 3.020 && m_jpsi    < 3.135"
# cuts_ += "&& minann_K  > 0.2"
# cuts_ += "&& ((psi_l0tos & 2) == 2)"
# cuts_ += " && MIPCHI2DV_k1 > 9 && MIPCHI2DV_k2 > 9 && MIPCHI2DV_k3 > 9"
# cuts_ += "&& ((psi_l1tos & 2) == 2)"
# cuts_ += "&& ((psi_l2tos & 2) == 2)"

# cuts_ += "&& ((psi3_l0tos & 2) == 2)"


# Weak cuts
# cuts_ = " DTFchi2ndof > 0"
# cuts_ += "&& DTFchi2ndof < 5"
# cuts_ += "&& DTFctau > 0.2"
# cuts_ += "&& pt_kaon[0] > 0.6 && pt_kaon[1] > 0.6" # && pt_kaon[2] > 0.6 "
# cuts_ += "&& m_jpsi    > 3.020 && m_jpsi    < 3.135"
# cuts_ += "&& minann_K  > 0.2"
# cuts_ += " && ((psi3_l0tos & 2) == 2)"


cuts_Bu = cuts_


def prntCuts(cuts, prefix=""):
    for cut in cuts.split("&&"):
        yield prefix + cut
