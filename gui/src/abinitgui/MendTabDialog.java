/*--
MendTabDialog.java - Created July 17, 2009

Copyright (c) 2009-2011 Flavio Miguel ABREU ARAUJO.
Université catholique de Louvain, Louvain-la-Neuve, Belgium
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

1. Redistributions of source code must retain the above copyright
notice, this list of conditions, and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright
notice, this list of conditions, and the disclaimer that follows
these conditions in the documentation and/or other materials
provided with the distribution.

3. The names of the author may not be used to endorse or promote
products derived from this software without specific prior written
permission.

In addition, we request (but do not require) that you include in the
end-user documentation provided with the redistribution and/or in the
software itself an acknowledgement equivalent to the following:
"This product includes software developed by the
Abinit Project (http://www.abinit.org/)."

THIS SOFTWARE IS PROVIDED ``AS IS'' AND ANY EXPRESSED OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED.  IN NO EVENT SHALL THE JDOM AUTHORS OR THE PROJECT
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
SUCH DAMAGE.

For more information on the Abinit Project, please see
<http://www.abinit.org/>.
 */
package abinitgui;

import java.awt.event.ActionListener;
import java.util.Enumeration;
import javax.swing.JCheckBox;
import javax.swing.JDialog;
import javax.swing.SwingUtilities;

//@SuppressWarnings("serial")
public class MendTabDialog extends JDialog {

    /** Creates new form MendTabDialog */
    public MendTabDialog(MainFrame parent, boolean modal, ActionListener al) {
        super(parent, modal);
        SwingUtilities.updateComponentTreeUI(this);
        initComponents();
        jB_H.addActionListener(al);
        jB_Sr.addActionListener(al);
        jB_Cm.addActionListener(al);
        jB_Am.addActionListener(al);
        jB_Eu.addActionListener(al);
        jB_Np.addActionListener(al);
        jB_Tb.addActionListener(al);
        jB_Cf.addActionListener(al);
        jB_Dy.addActionListener(al);
        jB_Gd.addActionListener(al);
        jB_Bk.addActionListener(al);
        jB_Ca.addActionListener(al);
        jB_Mg.addActionListener(al);
        jB_Be.addActionListener(al);
        jB_Sc.addActionListener(al);
        jB_Y.addActionListener(al);
        jB_Rf.addActionListener(al);
        jB_Zr.addActionListener(al);
        jB_Li.addActionListener(al);
        jB_Hf.addActionListener(al);
        jB_Ti.addActionListener(al);
        jB_Db.addActionListener(al);
        jB_Nb.addActionListener(al);
        jB_Ta.addActionListener(al);
        jB_V.addActionListener(al);
        jB_Sg.addActionListener(al);
        jB_Mo.addActionListener(al);
        jB_W.addActionListener(al);
        jB_Cr.addActionListener(al);
        jB_Na.addActionListener(al);
        jB_Tc.addActionListener(al);
        jB_Re.addActionListener(al);
        jB_Mn.addActionListener(al);
        jB_Fe.addActionListener(al);
        jB_Os.addActionListener(al);
        jB_Ru.addActionListener(al);
        jB_Co.addActionListener(al);
        jB_Ir.addActionListener(al);
        jB_Rh.addActionListener(al);
        jB_Ni.addActionListener(al);
        jB_K.addActionListener(al);
        jB_Pt.addActionListener(al);
        jB_Pd.addActionListener(al);
        jB_Cu.addActionListener(al);
        jB_Au.addActionListener(al);
        jB_Ag.addActionListener(al);
        jB_Zn.addActionListener(al);
        jB_Hg.addActionListener(al);
        jB_Cd.addActionListener(al);
        jB_B.addActionListener(al);
        jB_Ga.addActionListener(al);
        jB_Rb.addActionListener(al);
        jB_Al.addActionListener(al);
        jB_In.addActionListener(al);
        jB_Tl.addActionListener(al);
        jB_Sn.addActionListener(al);
        jB_Pb.addActionListener(al);
        jB_C.addActionListener(al);
        jB_Ge.addActionListener(al);
        jB_Si.addActionListener(al);
        jB_P.addActionListener(al);
        jB_As.addActionListener(al);
        jB_Cs.addActionListener(al);
        jB_N.addActionListener(al);
        jB_Sb.addActionListener(al);
        jB_Bi.addActionListener(al);
        jB_O.addActionListener(al);
        jB_Se.addActionListener(al);
        jB_S.addActionListener(al);
        jB_Po.addActionListener(al);
        jB_Te.addActionListener(al);
        jB_I.addActionListener(al);
        jB_At.addActionListener(al);
        jB_Fr.addActionListener(al);
        jB_F.addActionListener(al);
        jB_Br.addActionListener(al);
        jB_Cl.addActionListener(al);
        jB_La.addActionListener(al);
        jB_Rn.addActionListener(al);
        jB_Xe.addActionListener(al);
        jB_Kr.addActionListener(al);
        jB_Ar.addActionListener(al);
        jB_Ne.addActionListener(al);
        jB_He.addActionListener(al);
        jB_Ra.addActionListener(al);
        jB_Ac.addActionListener(al);
        jB_Ce.addActionListener(al);
        jB_Th.addActionListener(al);
        jB_Er.addActionListener(al);
        jB_Fm.addActionListener(al);
        jB_Ho.addActionListener(al);
        jB_Es.addActionListener(al);
        jB_Nd.addActionListener(al);
        jB_U.addActionListener(al);
        jB_Pr.addActionListener(al);
        jB_Ba.addActionListener(al);
        jB_Pa.addActionListener(al);
        jB_Yb.addActionListener(al);
        jB_Pm.addActionListener(al);
        jB_No.addActionListener(al);
        jB_Lr.addActionListener(al);
        jB_Lu.addActionListener(al);
        jB_Md.addActionListener(al);
        jB_Tm.addActionListener(al);
        jB_Sm.addActionListener(al);
        jB_Pu.addActionListener(al);

        jB_Bh.addActionListener(al);
        jB_Hs.addActionListener(al);

        jB_Uus.addActionListener(al);
        //jB_Uus.setVisible(false);
        jB_Uuo.addActionListener(al);
        //jB_Uuo.setVisible(false);

        jB_Mt.addActionListener(al);

        jB_Ds.addActionListener(al);
        //jB_Ds.setVisible(false);
        jB_Rg.addActionListener(al);
        //jB_Rg.setVisible(false);
        jB_Cn.addActionListener(al);
        //jB_Cn.setVisible(false);
        jB_Uut.addActionListener(al);
        //jB_Uut.setVisible(false);
        jB_Uuq.addActionListener(al);
        //jB_Uuq.setVisible(false);
        jB_Uup.addActionListener(al);
        //jB_Uup.setVisible(false);
        jB_Uuh.addActionListener(al);
        //jB_Uuh.setVisible(false);

        userPSPButton.addActionListener(al);
        GGA_FHI_CheckBox.addActionListener(al);
        GGA_HGH_CheckBox.addActionListener(al);
        LDA_Core_CheckBox.addActionListener(al);
        LDA_FHI_CheckBox.addActionListener(al);
        LDA_GTH_CheckBox.addActionListener(al);
        LDA_HGH_CheckBox.addActionListener(al);
        LDA_TM_CheckBox.addActionListener(al);
        LDA_Teter_CheckBox.addActionListener(al);
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    //@SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        buttonGroup = new javax.swing.ButtonGroup();
        jB_H = new javax.swing.JButton();
        jB_Li = new javax.swing.JButton();
        jB_Na = new javax.swing.JButton();
        jB_K = new javax.swing.JButton();
        jB_Rb = new javax.swing.JButton();
        jB_Cs = new javax.swing.JButton();
        jB_Fr = new javax.swing.JButton();
        jB_Ra = new javax.swing.JButton();
        jB_Ba = new javax.swing.JButton();
        jB_Sr = new javax.swing.JButton();
        jB_Ca = new javax.swing.JButton();
        jB_Mg = new javax.swing.JButton();
        jB_Be = new javax.swing.JButton();
        jB_Sc = new javax.swing.JButton();
        jB_Y = new javax.swing.JButton();
        jB_Rf = new javax.swing.JButton();
        jB_Zr = new javax.swing.JButton();
        jB_Hf = new javax.swing.JButton();
        jB_Ti = new javax.swing.JButton();
        jB_Db = new javax.swing.JButton();
        jB_Nb = new javax.swing.JButton();
        jB_Ta = new javax.swing.JButton();
        jB_V = new javax.swing.JButton();
        jB_Sg = new javax.swing.JButton();
        jB_Mo = new javax.swing.JButton();
        jB_W = new javax.swing.JButton();
        jB_Cr = new javax.swing.JButton();
        jB_Tc = new javax.swing.JButton();
        jB_Re = new javax.swing.JButton();
        jB_Mn = new javax.swing.JButton();
        jB_Fe = new javax.swing.JButton();
        jB_Os = new javax.swing.JButton();
        jB_Ru = new javax.swing.JButton();
        jB_Co = new javax.swing.JButton();
        jB_Ir = new javax.swing.JButton();
        jB_Rh = new javax.swing.JButton();
        jB_Ni = new javax.swing.JButton();
        jB_Pt = new javax.swing.JButton();
        jB_Pd = new javax.swing.JButton();
        jB_Cu = new javax.swing.JButton();
        jB_Au = new javax.swing.JButton();
        jB_Ag = new javax.swing.JButton();
        jB_Zn = new javax.swing.JButton();
        jB_Hg = new javax.swing.JButton();
        jB_Cd = new javax.swing.JButton();
        jB_B = new javax.swing.JButton();
        jB_Ga = new javax.swing.JButton();
        jB_Al = new javax.swing.JButton();
        jB_In = new javax.swing.JButton();
        jB_Tl = new javax.swing.JButton();
        jB_Sn = new javax.swing.JButton();
        jB_Pb = new javax.swing.JButton();
        jB_C = new javax.swing.JButton();
        jB_Ge = new javax.swing.JButton();
        jB_Si = new javax.swing.JButton();
        jB_P = new javax.swing.JButton();
        jB_As = new javax.swing.JButton();
        jB_N = new javax.swing.JButton();
        jB_Sb = new javax.swing.JButton();
        jB_Bi = new javax.swing.JButton();
        jB_O = new javax.swing.JButton();
        jB_Se = new javax.swing.JButton();
        jB_S = new javax.swing.JButton();
        jB_Po = new javax.swing.JButton();
        jB_Te = new javax.swing.JButton();
        jB_I = new javax.swing.JButton();
        jB_At = new javax.swing.JButton();
        jB_F = new javax.swing.JButton();
        jB_Br = new javax.swing.JButton();
        jB_Cl = new javax.swing.JButton();
        jB_Rn = new javax.swing.JButton();
        jB_Xe = new javax.swing.JButton();
        jB_Kr = new javax.swing.JButton();
        jB_Ar = new javax.swing.JButton();
        jB_Ne = new javax.swing.JButton();
        jB_He = new javax.swing.JButton();
        jB_La = new javax.swing.JButton();
        jB_Ac = new javax.swing.JButton();
        jB_Ce = new javax.swing.JButton();
        jB_Th = new javax.swing.JButton();
        jB_Er = new javax.swing.JButton();
        jB_Fm = new javax.swing.JButton();
        jB_Ho = new javax.swing.JButton();
        jB_Es = new javax.swing.JButton();
        jB_Nd = new javax.swing.JButton();
        jB_U = new javax.swing.JButton();
        jB_Pr = new javax.swing.JButton();
        jB_Pa = new javax.swing.JButton();
        jB_Yb = new javax.swing.JButton();
        jB_Pm = new javax.swing.JButton();
        jB_No = new javax.swing.JButton();
        jB_Lr = new javax.swing.JButton();
        jB_Lu = new javax.swing.JButton();
        jB_Md = new javax.swing.JButton();
        jB_Tm = new javax.swing.JButton();
        jB_Sm = new javax.swing.JButton();
        jB_Pu = new javax.swing.JButton();
        jB_Cm = new javax.swing.JButton();
        jB_Am = new javax.swing.JButton();
        jB_Eu = new javax.swing.JButton();
        jB_Np = new javax.swing.JButton();
        jB_Tb = new javax.swing.JButton();
        jB_Cf = new javax.swing.JButton();
        jB_Dy = new javax.swing.JButton();
        jB_Gd = new javax.swing.JButton();
        jB_Bk = new javax.swing.JButton();
        pspTypePanel = new javax.swing.JPanel();
        LDA_FHI_CheckBox = new javax.swing.JCheckBox();
        LDA_Core_CheckBox = new javax.swing.JCheckBox();
        LDA_TM_CheckBox = new javax.swing.JCheckBox();
        LDA_Teter_CheckBox = new javax.swing.JCheckBox();
        LDA_HGH_CheckBox = new javax.swing.JCheckBox();
        GGA_FHI_CheckBox = new javax.swing.JCheckBox();
        GGA_HGH_CheckBox = new javax.swing.JCheckBox();
        LDA_GTH_CheckBox = new javax.swing.JCheckBox();
        userPSPButton = new javax.swing.JButton();
        jB_Bh = new javax.swing.JButton();
        jB_Hs = new javax.swing.JButton();
        jB_Mt = new javax.swing.JButton();
        jB_Ds = new javax.swing.JButton();
        jB_Rg = new javax.swing.JButton();
        jB_Cn = new javax.swing.JButton();
        jB_Uut = new javax.swing.JButton();
        jB_Uuq = new javax.swing.JButton();
        jB_Uup = new javax.swing.JButton();
        jB_Uuh = new javax.swing.JButton();
        jB_Uus = new javax.swing.JButton();
        jB_Uuo = new javax.swing.JButton();

        setDefaultCloseOperation(javax.swing.WindowConstants.DISPOSE_ON_CLOSE);

        jB_H.setText("H");

        jB_Li.setText("Li");

        jB_Na.setText("Na");

        jB_K.setText("K");

        jB_Rb.setText("Rb");

        jB_Cs.setText("Cs");

        jB_Fr.setText("Fr");

        jB_Ra.setText("Ra");

        jB_Ba.setText("Ba");

        jB_Sr.setText("Sr");

        jB_Ca.setText("Ca");

        jB_Mg.setText("Mg");

        jB_Be.setText("Be");

        jB_Sc.setText("Sc");

        jB_Y.setText("Y");

        jB_Rf.setText("Rf");

        jB_Zr.setText("Zr");

        jB_Hf.setText("Hf");

        jB_Ti.setText("Ti");

        jB_Db.setText("Db");

        jB_Nb.setText("Nb");

        jB_Ta.setText("Ta");

        jB_V.setText("V");

        jB_Sg.setText("Sg");

        jB_Mo.setText("Mo");

        jB_W.setText("W");

        jB_Cr.setText("Cr");

        jB_Tc.setText("Tc");

        jB_Re.setText("Re");

        jB_Mn.setText("Mn");

        jB_Fe.setText("Fe");

        jB_Os.setText("Os");

        jB_Ru.setText("Ru");

        jB_Co.setText("Co");

        jB_Ir.setText("Ir");

        jB_Rh.setText("Rh");

        jB_Ni.setText("Ni");

        jB_Pt.setText("Pt");

        jB_Pd.setText("Pd");

        jB_Cu.setText("Cu");

        jB_Au.setText("Au");

        jB_Ag.setText("Ag");

        jB_Zn.setText("Zn");

        jB_Hg.setText("Hg");

        jB_Cd.setText("Cd");

        jB_B.setText("B");

        jB_Ga.setText("Ga");

        jB_Al.setText("Al");

        jB_In.setText("In");

        jB_Tl.setText("Tl");

        jB_Sn.setText("Sn");

        jB_Pb.setText("Pb");

        jB_C.setText("C");

        jB_Ge.setText("Ge");

        jB_Si.setText("Si");

        jB_P.setText("P");

        jB_As.setText("As");

        jB_N.setText("N");

        jB_Sb.setText("Sb");

        jB_Bi.setText("Bi");

        jB_O.setText("O");

        jB_Se.setText("Se");

        jB_S.setText("S");

        jB_Po.setText("Po");

        jB_Te.setText("Te");

        jB_I.setText("I");

        jB_At.setText("At");

        jB_F.setText("F");

        jB_Br.setText("Br");

        jB_Cl.setText("Cl");

        jB_Rn.setText("Rn");

        jB_Xe.setText("Xe");

        jB_Kr.setText("Kr");

        jB_Ar.setText("Ar");

        jB_Ne.setText("Ne");

        jB_He.setText("He");

        jB_La.setText("La");

        jB_Ac.setText("Ac");

        jB_Ce.setText("Ce");

        jB_Th.setText("Th");

        jB_Er.setText("Er");

        jB_Fm.setText("Fm");

        jB_Ho.setText("Ho");

        jB_Es.setText("Es");

        jB_Nd.setText("Nd");

        jB_U.setText("U");

        jB_Pr.setText("Pr");

        jB_Pa.setText("Pa");

        jB_Yb.setText("Yb");

        jB_Pm.setText("Pm");

        jB_No.setText("No");

        jB_Lr.setText("Lr");

        jB_Lu.setText("Lu");

        jB_Md.setText("Md");

        jB_Tm.setText("Tm");

        jB_Sm.setText("Sm");

        jB_Pu.setText("Pu");

        jB_Cm.setText("Cm");

        jB_Am.setText("Am");

        jB_Eu.setText("Eu");

        jB_Np.setText("Np");

        jB_Tb.setText("Tb");

        jB_Cf.setText("Cf");

        jB_Dy.setText("Dy");

        jB_Gd.setText("Gd");

        jB_Bk.setText("Bk");

        pspTypePanel.setBorder(javax.swing.BorderFactory.createTitledBorder(null, "Choose pseudopotential type ...", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("Arial", 3, 14), java.awt.Color.blue)); // NOI18N

        buttonGroup.add(LDA_FHI_CheckBox);
        LDA_FHI_CheckBox.setText("LDA FHI");
        LDA_FHI_CheckBox.setActionCommand("LDA_FHI");
        LDA_FHI_CheckBox.setName("LDA_FHI"); // NOI18N

        buttonGroup.add(LDA_Core_CheckBox);
        LDA_Core_CheckBox.setText("LDA Core");
        LDA_Core_CheckBox.setActionCommand("LDA_Core");
        LDA_Core_CheckBox.setName("LDA_Core"); // NOI18N

        buttonGroup.add(LDA_TM_CheckBox);
        LDA_TM_CheckBox.setSelected(true);
        LDA_TM_CheckBox.setText("LDA TM");
        LDA_TM_CheckBox.setActionCommand("LDA_TM");
        LDA_TM_CheckBox.setName("LDA_TM"); // NOI18N
        LDA_TM_CheckBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                LDA_TM_CheckBoxActionPerformed(evt);
            }
        });

        buttonGroup.add(LDA_Teter_CheckBox);
        LDA_Teter_CheckBox.setText("LDA Teter");
        LDA_Teter_CheckBox.setActionCommand("LDA_Teter");
        LDA_Teter_CheckBox.setEnabled(false);
        LDA_Teter_CheckBox.setName("LDA_Teter"); // NOI18N

        buttonGroup.add(LDA_HGH_CheckBox);
        LDA_HGH_CheckBox.setText("LDA HGH");
        LDA_HGH_CheckBox.setActionCommand("LDA_HGH");
        LDA_HGH_CheckBox.setEnabled(false);
        LDA_HGH_CheckBox.setName("LDA_HGH"); // NOI18N

        buttonGroup.add(GGA_FHI_CheckBox);
        GGA_FHI_CheckBox.setText("GGA FHI");
        GGA_FHI_CheckBox.setActionCommand("GGA_FHI");
        GGA_FHI_CheckBox.setName("GGA_FHI"); // NOI18N

        buttonGroup.add(GGA_HGH_CheckBox);
        GGA_HGH_CheckBox.setText("GGA HGH");
        GGA_HGH_CheckBox.setActionCommand("GGA_HGH");
        GGA_HGH_CheckBox.setName("GGA_HGH"); // NOI18N

        buttonGroup.add(LDA_GTH_CheckBox);
        LDA_GTH_CheckBox.setText("LDA GTH");
        LDA_GTH_CheckBox.setActionCommand("LDA_GTH");
        LDA_GTH_CheckBox.setName("LDA_GTH"); // NOI18N

        javax.swing.GroupLayout pspTypePanelLayout = new javax.swing.GroupLayout(pspTypePanel);
        pspTypePanel.setLayout(pspTypePanelLayout);
        pspTypePanelLayout.setHorizontalGroup(
            pspTypePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pspTypePanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pspTypePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(LDA_TM_CheckBox)
                    .addComponent(LDA_GTH_CheckBox))
                .addGap(18, 18, 18)
                .addGroup(pspTypePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(LDA_FHI_CheckBox)
                    .addComponent(GGA_FHI_CheckBox))
                .addGap(18, 18, 18)
                .addGroup(pspTypePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(GGA_HGH_CheckBox)
                    .addComponent(LDA_Core_CheckBox))
                .addGap(18, 18, 18)
                .addGroup(pspTypePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addComponent(LDA_HGH_CheckBox)
                    .addComponent(LDA_Teter_CheckBox))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        pspTypePanelLayout.setVerticalGroup(
            pspTypePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(pspTypePanelLayout.createSequentialGroup()
                .addContainerGap()
                .addGroup(pspTypePanelLayout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(pspTypePanelLayout.createSequentialGroup()
                        .addComponent(LDA_FHI_CheckBox)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(GGA_FHI_CheckBox))
                    .addGroup(pspTypePanelLayout.createSequentialGroup()
                        .addComponent(LDA_TM_CheckBox)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(LDA_GTH_CheckBox))
                    .addGroup(pspTypePanelLayout.createSequentialGroup()
                        .addComponent(LDA_Core_CheckBox)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(GGA_HGH_CheckBox))
                    .addGroup(pspTypePanelLayout.createSequentialGroup()
                        .addComponent(LDA_Teter_CheckBox)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(LDA_HGH_CheckBox)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        LDA_FHI_CheckBox.getAccessibleContext().setAccessibleName("LDA_FHI");
        LDA_Core_CheckBox.getAccessibleContext().setAccessibleName("LDA_Core");
        LDA_TM_CheckBox.getAccessibleContext().setAccessibleName("LDA_TM");
        LDA_Teter_CheckBox.getAccessibleContext().setAccessibleName("LDA_Teter");
        LDA_HGH_CheckBox.getAccessibleContext().setAccessibleName("LDA_HGH");
        GGA_FHI_CheckBox.getAccessibleContext().setAccessibleName("GGA_FHI");
        GGA_HGH_CheckBox.getAccessibleContext().setAccessibleName("GGA_HGH");
        LDA_GTH_CheckBox.getAccessibleContext().setAccessibleName("LDA_GTH");

        userPSPButton.setText("User PSP");
        userPSPButton.setActionCommand("UserPSP");

        jB_Bh.setText("Bh");

        jB_Hs.setText("Hs");

        jB_Mt.setText("Mt");

        jB_Ds.setText("Ds");

        jB_Rg.setText("Rg");

        jB_Cn.setText("Cn");

        jB_Uut.setText("Uut");

        jB_Uuq.setText("Uuq");

        jB_Uup.setText("Uup");

        jB_Uuh.setText("Uuh");

        jB_Uus.setText("Uus");

        jB_Uuo.setText("Uuo");

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(pspTypePanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(344, 344, 344)
                                .addComponent(userPSPButton))
                            .addGroup(layout.createSequentialGroup()
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(layout.createSequentialGroup()
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(jB_Li)
                                            .addComponent(jB_Na)
                                            .addComponent(jB_K)
                                            .addComponent(jB_Rb)
                                            .addComponent(jB_Cs)
                                            .addComponent(jB_Fr))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(jB_Be)
                                            .addComponent(jB_Mg)
                                            .addGroup(layout.createSequentialGroup()
                                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                    .addComponent(jB_Ca)
                                                    .addComponent(jB_Sr)
                                                    .addComponent(jB_Ba)
                                                    .addComponent(jB_Ra))
                                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                    .addComponent(jB_Sc)
                                                    .addComponent(jB_Y))
                                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                    .addComponent(jB_Ti)
                                                    .addComponent(jB_Zr)
                                                    .addComponent(jB_Hf)
                                                    .addComponent(jB_Rf))
                                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                    .addComponent(jB_V)
                                                    .addComponent(jB_Nb)
                                                    .addComponent(jB_Ta)
                                                    .addComponent(jB_Db))
                                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                    .addGroup(layout.createSequentialGroup()
                                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                            .addComponent(jB_Cr)
                                                            .addComponent(jB_Mo)
                                                            .addComponent(jB_W))
                                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                            .addComponent(jB_Mn)
                                                            .addComponent(jB_Tc)
                                                            .addComponent(jB_Re))
                                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                            .addComponent(jB_Fe)
                                                            .addComponent(jB_Ru)
                                                            .addComponent(jB_Os))
                                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                            .addComponent(jB_Co)
                                                            .addComponent(jB_Rh)
                                                            .addComponent(jB_Ir))
                                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                            .addComponent(jB_Ni)
                                                            .addComponent(jB_Pd)
                                                            .addComponent(jB_Pt))
                                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                            .addComponent(jB_Cu)
                                                            .addComponent(jB_Ag)
                                                            .addComponent(jB_Au))
                                                        .addGap(6, 6, 6)
                                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                                            .addComponent(jB_Zn)
                                                            .addComponent(jB_Cd)
                                                            .addComponent(jB_Hg)))
                                                    .addGroup(layout.createSequentialGroup()
                                                        .addComponent(jB_Sg)
                                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                        .addComponent(jB_Bh)
                                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                        .addComponent(jB_Hs)
                                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                        .addComponent(jB_Mt)
                                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                        .addComponent(jB_Ds)
                                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                        .addComponent(jB_Rg)
                                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                        .addComponent(jB_Cn))))))
                                    .addComponent(jB_H))
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(layout.createSequentialGroup()
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(jB_B)
                                            .addComponent(jB_Al)
                                            .addComponent(jB_Ga)
                                            .addComponent(jB_In)
                                            .addComponent(jB_Tl))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(jB_C)
                                            .addComponent(jB_Si)
                                            .addComponent(jB_Ge)
                                            .addComponent(jB_Sn)
                                            .addComponent(jB_Pb))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(jB_N)
                                            .addComponent(jB_P)
                                            .addComponent(jB_As)
                                            .addComponent(jB_Sb)
                                            .addComponent(jB_Bi))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(jB_O)
                                            .addComponent(jB_S)
                                            .addComponent(jB_Se)
                                            .addComponent(jB_Te)
                                            .addComponent(jB_Po))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(jB_F)
                                            .addComponent(jB_Cl)
                                            .addComponent(jB_Br)
                                            .addComponent(jB_I)
                                            .addComponent(jB_At))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addComponent(jB_Ne)
                                            .addComponent(jB_Ar)
                                            .addComponent(jB_Kr)
                                            .addComponent(jB_Xe)
                                            .addComponent(jB_Rn)
                                            .addComponent(jB_He)))
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(jB_Uut)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Uuq)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Uup)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Uuh)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Uus)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Uuo))))))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(136, 136, 136)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_La)
                            .addComponent(jB_Ac))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Ce)
                            .addComponent(jB_Th))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Pr)
                            .addComponent(jB_Pa))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Nd)
                            .addComponent(jB_U))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Pm)
                            .addComponent(jB_Np))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Sm)
                            .addComponent(jB_Pu))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Eu)
                            .addComponent(jB_Am))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Gd)
                            .addComponent(jB_Cm))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Tb)
                            .addComponent(jB_Bk))
                        .addGap(6, 6, 6)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Dy)
                            .addComponent(jB_Cf))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Ho)
                            .addComponent(jB_Es))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Er)
                            .addComponent(jB_Fm))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Tm)
                            .addComponent(jB_Md))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Yb)
                            .addComponent(jB_No))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addComponent(jB_Lu)
                            .addComponent(jB_Lr))))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        layout.linkSize(javax.swing.SwingConstants.HORIZONTAL, new java.awt.Component[] {jB_Ac, jB_Ag, jB_Al, jB_Am, jB_Ar, jB_As, jB_At, jB_Au, jB_B, jB_Ba, jB_Be, jB_Bh, jB_Bi, jB_Bk, jB_Br, jB_C, jB_Ca, jB_Cd, jB_Ce, jB_Cf, jB_Cl, jB_Cm, jB_Cn, jB_Co, jB_Cr, jB_Cs, jB_Cu, jB_Db, jB_Ds, jB_Dy, jB_Er, jB_Es, jB_Eu, jB_F, jB_Fe, jB_Fm, jB_Fr, jB_Ga, jB_Gd, jB_Ge, jB_H, jB_He, jB_Hf, jB_Hg, jB_Ho, jB_Hs, jB_I, jB_In, jB_Ir, jB_K, jB_Kr, jB_La, jB_Li, jB_Lr, jB_Lu, jB_Md, jB_Mg, jB_Mn, jB_Mo, jB_Mt, jB_N, jB_Na, jB_Nb, jB_Nd, jB_Ne, jB_Ni, jB_No, jB_Np, jB_O, jB_Os, jB_P, jB_Pa, jB_Pb, jB_Pd, jB_Pm, jB_Po, jB_Pr, jB_Pt, jB_Pu, jB_Ra, jB_Rb, jB_Re, jB_Rf, jB_Rg, jB_Rh, jB_Rn, jB_Ru, jB_S, jB_Sb, jB_Sc, jB_Se, jB_Sg, jB_Si, jB_Sm, jB_Sn, jB_Sr, jB_Ta, jB_Tb, jB_Tc, jB_Te, jB_Th, jB_Ti, jB_Tl, jB_Tm, jB_U, jB_Uuh, jB_Uuo, jB_Uup, jB_Uuq, jB_Uus, jB_Uut, jB_V, jB_W, jB_Xe, jB_Y, jB_Yb, jB_Zn, jB_Zr});

        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_He)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jB_C)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Si)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Ge)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Sn)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Pb))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jB_B)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Al)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Ga)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_In)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Tl))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jB_N)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_P)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_As)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Sb)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Bi))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jB_O)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_S)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Se)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Te)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Po))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jB_F)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Cl)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Br)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_I)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_At))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jB_Ne)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Ar)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Kr)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Xe)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Rn))))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_H)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jB_Li)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Na)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_K)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Rb)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Cs)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Fr))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jB_Be)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addComponent(jB_Mg)
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(jB_Sc)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Y))
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(jB_Ca)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Sr)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Ba)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Ra))
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(jB_Ti)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Zr)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Hf)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Rf))
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(jB_V)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Nb)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Ta)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Db))
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(jB_Fe)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Ru)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Os))
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(jB_Co)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Rh)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Ir))
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(jB_Ni)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Pd)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Pt))
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(jB_Cu)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Ag)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Au))
                                    .addGroup(layout.createSequentialGroup()
                                        .addComponent(jB_Zn)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Cd)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addComponent(jB_Hg))
                                    .addGroup(layout.createSequentialGroup()
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                                            .addGroup(layout.createSequentialGroup()
                                                .addComponent(jB_Cr)
                                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                .addComponent(jB_Mo)
                                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                .addComponent(jB_W))
                                            .addGroup(layout.createSequentialGroup()
                                                .addComponent(jB_Mn)
                                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                .addComponent(jB_Tc)
                                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                                .addComponent(jB_Re)))
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                                            .addComponent(jB_Sg)
                                            .addComponent(jB_Bh)
                                            .addComponent(jB_Hs)
                                            .addComponent(jB_Mt)
                                            .addComponent(jB_Ds)
                                            .addComponent(jB_Rg)
                                            .addComponent(jB_Cn)
                                            .addComponent(jB_Uut)
                                            .addComponent(jB_Uuq)
                                            .addComponent(jB_Uup)
                                            .addComponent(jB_Uuh)
                                            .addComponent(jB_Uus)
                                            .addComponent(jB_Uuo))))))))
                .addGap(18, 18, 18)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Er)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Fm))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Ho)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Es))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Tm)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Md))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Yb)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_No))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Lu)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Lr))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_La)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Ac))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Ce)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Th))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Pr)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Pa))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Sm)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Pu))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Eu)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Am))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Gd)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Cm))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Tb)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Bk))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Dy)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Cf))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Nd)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_U))
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(jB_Pm)
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(jB_Np)))
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(12, 12, 12)
                        .addComponent(pspTypePanel, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                    .addGroup(layout.createSequentialGroup()
                        .addGap(55, 55, 55)
                        .addComponent(userPSPButton)))
                .addContainerGap(javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        layout.linkSize(javax.swing.SwingConstants.VERTICAL, new java.awt.Component[] {jB_Ac, jB_Ag, jB_Al, jB_Am, jB_Ar, jB_As, jB_At, jB_Au, jB_B, jB_Ba, jB_Be, jB_Bh, jB_Bi, jB_Bk, jB_Br, jB_C, jB_Ca, jB_Cd, jB_Ce, jB_Cf, jB_Cl, jB_Cm, jB_Cn, jB_Co, jB_Cr, jB_Cs, jB_Cu, jB_Db, jB_Ds, jB_Dy, jB_Er, jB_Es, jB_Eu, jB_F, jB_Fe, jB_Fm, jB_Fr, jB_Ga, jB_Gd, jB_Ge, jB_H, jB_He, jB_Hf, jB_Hg, jB_Ho, jB_Hs, jB_I, jB_In, jB_Ir, jB_K, jB_Kr, jB_La, jB_Li, jB_Lr, jB_Lu, jB_Md, jB_Mg, jB_Mn, jB_Mo, jB_Mt, jB_N, jB_Na, jB_Nb, jB_Nd, jB_Ne, jB_Ni, jB_No, jB_Np, jB_O, jB_Os, jB_P, jB_Pa, jB_Pb, jB_Pd, jB_Pm, jB_Po, jB_Pr, jB_Pt, jB_Pu, jB_Ra, jB_Rb, jB_Re, jB_Rf, jB_Rg, jB_Rh, jB_Rn, jB_Ru, jB_S, jB_Sb, jB_Sc, jB_Se, jB_Sg, jB_Si, jB_Sm, jB_Sn, jB_Sr, jB_Ta, jB_Tb, jB_Tc, jB_Te, jB_Th, jB_Ti, jB_Tl, jB_Tm, jB_U, jB_Uuh, jB_Uuo, jB_Uup, jB_Uuq, jB_Uus, jB_Uut, jB_V, jB_W, jB_Xe, jB_Y, jB_Yb, jB_Zn, jB_Zr});

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void LDA_TM_CheckBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_LDA_TM_CheckBoxActionPerformed
        // TODO add your handling code here:
    }//GEN-LAST:event_LDA_TM_CheckBoxActionPerformed

    public String getPSPSelected() {
        Enumeration buttons = buttonGroup.getElements();
        while (buttons.hasMoreElements()) {
            JCheckBox cb = (JCheckBox) buttons.nextElement();
            if (cb.isSelected()) {
                return cb.getName();
            }
        }
        return null;
    }

    public javax.swing.JButton getMendButton(String atom) {
        if (atom.equals("H")) {
            return jB_H;
        }
        if (atom.equals("Sr")) {
            return jB_Sr;
        }
        if (atom.equals("Cm")) {
            return jB_Cm;
        }
        if (atom.equals("Am")) {
            return jB_Am;
        }
        if (atom.equals("Eu")) {
            return jB_Eu;
        }
        if (atom.equals("Np")) {
            return jB_Np;
        }
        if (atom.equals("Tb")) {
            return jB_Tb;
        }
        if (atom.equals("Cf")) {
            return jB_Cf;
        }
        if (atom.equals("Dy")) {
            return jB_Dy;
        }
        if (atom.equals("Gd")) {
            return jB_Gd;
        }
        if (atom.equals("Bk")) {
            return jB_Bk;
        }
        if (atom.equals("Ca")) {
            return jB_Ca;
        }
        if (atom.equals("Mg")) {
            return jB_Mg;
        }
        if (atom.equals("Be")) {
            return jB_Be;
        }
        if (atom.equals("Sc")) {
            return jB_Sc;
        }
        if (atom.equals("Y")) {
            return jB_Y;
        }
        if (atom.equals("Rf")) {
            return jB_Rf;
        }
        if (atom.equals("Zr")) {
            return jB_Zr;
        }
        if (atom.equals("Li")) {
            return jB_Li;
        }
        if (atom.equals("Hf")) {
            return jB_Hf;
        }
        if (atom.equals("Ti")) {
            return jB_Ti;
        }
        if (atom.equals("Db")) {
            return jB_Db;
        }
        if (atom.equals("Nb")) {
            return jB_Nb;
        }
        if (atom.equals("Ta")) {
            return jB_Ta;
        }
        if (atom.equals("V")) {
            return jB_V;
        }
        if (atom.equals("Sg")) {
            return jB_Sg;
        }
        if (atom.equals("Mo")) {
            return jB_Mo;
        }
        if (atom.equals("W")) {
            return jB_W;
        }
        if (atom.equals("Cr")) {
            return jB_Cr;
        }
        if (atom.equals("Na")) {
            return jB_Na;
        }
        if (atom.equals("Tc")) {
            return jB_Tc;
        }
        if (atom.equals("Re")) {
            return jB_Re;
        }
        if (atom.equals("Mn")) {
            return jB_Mn;
        }
        if (atom.equals("Fe")) {
            return jB_Fe;
        }
        if (atom.equals("Os")) {
            return jB_Os;
        }
        if (atom.equals("Ru")) {
            return jB_Ru;
        }
        if (atom.equals("Co")) {
            return jB_Co;
        }
        if (atom.equals("Ir")) {
            return jB_Ir;
        }
        if (atom.equals("Rh")) {
            return jB_Rh;
        }
        if (atom.equals("Ni")) {
            return jB_Ni;
        }
        if (atom.equals("K")) {
            return jB_K;
        }
        if (atom.equals("Pt")) {
            return jB_Pt;
        }
        if (atom.equals("Pd")) {
            return jB_Pd;
        }
        if (atom.equals("Cu")) {
            return jB_Cu;
        }
        if (atom.equals("Au")) {
            return jB_Au;
        }
        if (atom.equals("Ag")) {
            return jB_Ag;
        }
        if (atom.equals("Zn")) {
            return jB_Zn;
        }
        if (atom.equals("Hg")) {
            return jB_Hg;
        }
        if (atom.equals("Cd")) {
            return jB_Cd;
        }
        if (atom.equals("B")) {
            return jB_B;
        }
        if (atom.equals("Ga")) {
            return jB_Ga;
        }
        if (atom.equals("Rb")) {
            return jB_Rb;
        }
        if (atom.equals("Al")) {
            return jB_Al;
        }
        if (atom.equals("In")) {
            return jB_In;
        }
        if (atom.equals("Tl")) {
            return jB_Tl;
        }
        if (atom.equals("Sn")) {
            return jB_Sn;
        }
        if (atom.equals("Pb")) {
            return jB_Pb;
        }
        if (atom.equals("C")) {
            return jB_C;
        }
        if (atom.equals("Ge")) {
            return jB_Ge;
        }
        if (atom.equals("Si")) {
            return jB_Si;
        }
        if (atom.equals("P")) {
            return jB_P;
        }
        if (atom.equals("As")) {
            return jB_As;
        }
        if (atom.equals("Cs")) {
            return jB_Cs;
        }
        if (atom.equals("N")) {
            return jB_N;
        }
        if (atom.equals("Sb")) {
            return jB_Sb;
        }
        if (atom.equals("Bi")) {
            return jB_Bi;
        }
        if (atom.equals("O")) {
            return jB_O;
        }
        if (atom.equals("Se")) {
            return jB_Se;
        }
        if (atom.equals("S")) {
            return jB_S;
        }
        if (atom.equals("Po")) {
            return jB_Po;
        }
        if (atom.equals("Te")) {
            return jB_Te;
        }
        if (atom.equals("I")) {
            return jB_I;
        }
        if (atom.equals("At")) {
            return jB_At;
        }
        if (atom.equals("Fr")) {
            return jB_Fr;
        }
        if (atom.equals("F")) {
            return jB_F;
        }
        if (atom.equals("Br")) {
            return jB_Br;
        }
        if (atom.equals("Cl")) {
            return jB_Cl;
        }
        if (atom.equals("La")) {
            return jB_La;
        }
        if (atom.equals("Rn")) {
            return jB_Rn;
        }
        if (atom.equals("Xe")) {
            return jB_Xe;
        }
        if (atom.equals("Kr")) {
            return jB_Kr;
        }
        if (atom.equals("Ar")) {
            return jB_Ar;
        }
        if (atom.equals("Ne")) {
            return jB_Ne;
        }
        if (atom.equals("He")) {
            return jB_He;
        }
        if (atom.equals("Ra")) {
            return jB_Ra;
        }
        if (atom.equals("Ac")) {
            return jB_Ac;
        }
        if (atom.equals("Ce")) {
            return jB_Ce;
        }
        if (atom.equals("Th")) {
            return jB_Th;
        }
        if (atom.equals("Er")) {
            return jB_Er;
        }
        if (atom.equals("Fm")) {
            return jB_Fm;
        }
        if (atom.equals("Ho")) {
            return jB_Ho;
        }
        if (atom.equals("Es")) {
            return jB_Es;
        }
        if (atom.equals("Nd")) {
            return jB_Nd;
        }
        if (atom.equals("U")) {
            return jB_U;
        }
        if (atom.equals("Pr")) {
            return jB_Pr;
        }
        if (atom.equals("Ba")) {
            return jB_Ba;
        }
        if (atom.equals("Pa")) {
            return jB_Pa;
        }
        if (atom.equals("Yb")) {
            return jB_Yb;
        }
        if (atom.equals("Pm")) {
            return jB_Pm;
        }
        if (atom.equals("No")) {
            return jB_No;
        }
        if (atom.equals("Lr")) {
            return jB_Lr;
        }
        if (atom.equals("Lu")) {
            return jB_Lu;
        }
        if (atom.equals("Md")) {
            return jB_Md;
        }
        if (atom.equals("Tm")) {
            return jB_Tm;
        }
        if (atom.equals("Sm")) {
            return jB_Sm;
        }
        if (atom.equals("Pu")) {
            return jB_Pu;
        }
        if (atom.equals("Bh")) {
            return jB_Bh;
        }
        if (atom.equals("Hs")) {
            return jB_Hs;
        }
        if (atom.equals("Uus")) {
            return jB_Uus;
        }
        if (atom.equals("Uuo")) {
            return jB_Uuo;
        }
        if (atom.equals("Mt")) {
            return jB_Mt;
        }
        if (atom.equals("Ds")) {
            return jB_Ds;
        }
        if (atom.equals("Rg")) {
            return jB_Rg;
        }
        if (atom.equals("Cn")) {
            return jB_Cn;
        }
        if (atom.equals("Uut")) {
            return jB_Uut;
        }
        if (atom.equals("Uuq")) {
            return jB_Uuq;
        }
        if (atom.equals("Uup")) {
            return jB_Uup;
        }
        if (atom.equals("Uuh")) {
            return jB_Uuh;
        } else {
            return null;
        }
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JCheckBox GGA_FHI_CheckBox;
    private javax.swing.JCheckBox GGA_HGH_CheckBox;
    private javax.swing.JCheckBox LDA_Core_CheckBox;
    private javax.swing.JCheckBox LDA_FHI_CheckBox;
    private javax.swing.JCheckBox LDA_GTH_CheckBox;
    private javax.swing.JCheckBox LDA_HGH_CheckBox;
    private javax.swing.JCheckBox LDA_TM_CheckBox;
    private javax.swing.JCheckBox LDA_Teter_CheckBox;
    private javax.swing.ButtonGroup buttonGroup;
    private javax.swing.JButton jB_Ac;
    private javax.swing.JButton jB_Ag;
    private javax.swing.JButton jB_Al;
    private javax.swing.JButton jB_Am;
    private javax.swing.JButton jB_Ar;
    private javax.swing.JButton jB_As;
    private javax.swing.JButton jB_At;
    private javax.swing.JButton jB_Au;
    private javax.swing.JButton jB_B;
    private javax.swing.JButton jB_Ba;
    private javax.swing.JButton jB_Be;
    private javax.swing.JButton jB_Bh;
    private javax.swing.JButton jB_Bi;
    private javax.swing.JButton jB_Bk;
    private javax.swing.JButton jB_Br;
    private javax.swing.JButton jB_C;
    private javax.swing.JButton jB_Ca;
    private javax.swing.JButton jB_Cd;
    private javax.swing.JButton jB_Ce;
    private javax.swing.JButton jB_Cf;
    private javax.swing.JButton jB_Cl;
    private javax.swing.JButton jB_Cm;
    private javax.swing.JButton jB_Cn;
    private javax.swing.JButton jB_Co;
    private javax.swing.JButton jB_Cr;
    private javax.swing.JButton jB_Cs;
    private javax.swing.JButton jB_Cu;
    private javax.swing.JButton jB_Db;
    private javax.swing.JButton jB_Ds;
    private javax.swing.JButton jB_Dy;
    private javax.swing.JButton jB_Er;
    private javax.swing.JButton jB_Es;
    private javax.swing.JButton jB_Eu;
    private javax.swing.JButton jB_F;
    private javax.swing.JButton jB_Fe;
    private javax.swing.JButton jB_Fm;
    private javax.swing.JButton jB_Fr;
    private javax.swing.JButton jB_Ga;
    private javax.swing.JButton jB_Gd;
    private javax.swing.JButton jB_Ge;
    private javax.swing.JButton jB_H;
    private javax.swing.JButton jB_He;
    private javax.swing.JButton jB_Hf;
    private javax.swing.JButton jB_Hg;
    private javax.swing.JButton jB_Ho;
    private javax.swing.JButton jB_Hs;
    private javax.swing.JButton jB_I;
    private javax.swing.JButton jB_In;
    private javax.swing.JButton jB_Ir;
    private javax.swing.JButton jB_K;
    private javax.swing.JButton jB_Kr;
    private javax.swing.JButton jB_La;
    private javax.swing.JButton jB_Li;
    private javax.swing.JButton jB_Lr;
    private javax.swing.JButton jB_Lu;
    private javax.swing.JButton jB_Md;
    private javax.swing.JButton jB_Mg;
    private javax.swing.JButton jB_Mn;
    private javax.swing.JButton jB_Mo;
    private javax.swing.JButton jB_Mt;
    private javax.swing.JButton jB_N;
    private javax.swing.JButton jB_Na;
    private javax.swing.JButton jB_Nb;
    private javax.swing.JButton jB_Nd;
    private javax.swing.JButton jB_Ne;
    private javax.swing.JButton jB_Ni;
    private javax.swing.JButton jB_No;
    private javax.swing.JButton jB_Np;
    private javax.swing.JButton jB_O;
    private javax.swing.JButton jB_Os;
    private javax.swing.JButton jB_P;
    private javax.swing.JButton jB_Pa;
    private javax.swing.JButton jB_Pb;
    private javax.swing.JButton jB_Pd;
    private javax.swing.JButton jB_Pm;
    private javax.swing.JButton jB_Po;
    private javax.swing.JButton jB_Pr;
    private javax.swing.JButton jB_Pt;
    private javax.swing.JButton jB_Pu;
    private javax.swing.JButton jB_Ra;
    private javax.swing.JButton jB_Rb;
    private javax.swing.JButton jB_Re;
    private javax.swing.JButton jB_Rf;
    private javax.swing.JButton jB_Rg;
    private javax.swing.JButton jB_Rh;
    private javax.swing.JButton jB_Rn;
    private javax.swing.JButton jB_Ru;
    private javax.swing.JButton jB_S;
    private javax.swing.JButton jB_Sb;
    private javax.swing.JButton jB_Sc;
    private javax.swing.JButton jB_Se;
    private javax.swing.JButton jB_Sg;
    private javax.swing.JButton jB_Si;
    private javax.swing.JButton jB_Sm;
    private javax.swing.JButton jB_Sn;
    private javax.swing.JButton jB_Sr;
    private javax.swing.JButton jB_Ta;
    private javax.swing.JButton jB_Tb;
    private javax.swing.JButton jB_Tc;
    private javax.swing.JButton jB_Te;
    private javax.swing.JButton jB_Th;
    private javax.swing.JButton jB_Ti;
    private javax.swing.JButton jB_Tl;
    private javax.swing.JButton jB_Tm;
    private javax.swing.JButton jB_U;
    private javax.swing.JButton jB_Uuh;
    private javax.swing.JButton jB_Uuo;
    private javax.swing.JButton jB_Uup;
    private javax.swing.JButton jB_Uuq;
    private javax.swing.JButton jB_Uus;
    private javax.swing.JButton jB_Uut;
    private javax.swing.JButton jB_V;
    private javax.swing.JButton jB_W;
    private javax.swing.JButton jB_Xe;
    private javax.swing.JButton jB_Y;
    private javax.swing.JButton jB_Yb;
    private javax.swing.JButton jB_Zn;
    private javax.swing.JButton jB_Zr;
    private javax.swing.JPanel pspTypePanel;
    private javax.swing.JButton userPSPButton;
    // End of variables declaration//GEN-END:variables
}
