/*--
 AboutDialog.java - Created September 2, 2009

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

//@SuppressWarnings("serial")
public class AboutDialog extends javax.swing.JDialog {

    /** Creates new form AboutDialog */
    public AboutDialog(java.awt.Frame parent, boolean modal) {
        super(parent, modal);
        initComponents();
        versionLabel.setText("Version: 0.1.0.4 (April 2011)");
        osLabel.setText("Operating System: " + Utils.osName() + " (" + Utils.osArch() + ")");
        javaVersionLabel.setText("Java version number: " + Utils.javaVersion());
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    @SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        jPanel1 = new javax.swing.JPanel();
        abinitBannerLabel = new javax.swing.JLabel();
        versionLabel = new javax.swing.JLabel();
        authorLabel = new javax.swing.JLabel();
        emailLabel = new javax.swing.JLabel();
        osLabel = new javax.swing.JLabel();
        abinitVersionLabel = new javax.swing.JLabel();
        javaVersionLabel = new javax.swing.JLabel();
        emailLabel1 = new javax.swing.JLabel();

        setTitle("About");
        setResizable(false);

        jPanel1.setBackground(java.awt.Color.white);

        abinitBannerLabel.setFont(new java.awt.Font("Arial", 1, 70));
        abinitBannerLabel.setText("AbinitGUI");

        versionLabel.setFont(new java.awt.Font("Serif", 1, 14));
        versionLabel.setText("Version: ??");

        authorLabel.setFont(new java.awt.Font("Serif", 1, 14));
        authorLabel.setText("Developed by: Flavio ABREU ARAUJO");

        emailLabel.setFont(new java.awt.Font("Serif", 1, 14));
        emailLabel.setText("email: abinitgui@gmail.com");

        osLabel.setFont(new java.awt.Font("Serif", 1, 14));
        osLabel.setText("Operating System: ??");

        abinitVersionLabel.setFont(new java.awt.Font("Serif", 1, 14));
        abinitVersionLabel.setText("Designed for ABINIT version: 6+");

        javaVersionLabel.setFont(new java.awt.Font("Serif", 1, 14));
        javaVersionLabel.setText("Java version number: ??");

        emailLabel1.setFont(new java.awt.Font("Serif", 1, 14)); // NOI18N
        emailLabel1.setForeground(new java.awt.Color(246, 10, 13));
        emailLabel1.setText("<HTML><b>This product comes without any warranty. The potential damages are entirely the responsibility of the user.</b></HTML>");
        emailLabel1.setVerticalAlignment(javax.swing.SwingConstants.TOP);
        emailLabel1.setFocusable(false);
        emailLabel1.setHorizontalTextPosition(javax.swing.SwingConstants.LEFT);
        emailLabel1.setInheritsPopupMenu(false);
        emailLabel1.setVerticalTextPosition(javax.swing.SwingConstants.TOP);

        javax.swing.GroupLayout jPanel1Layout = new javax.swing.GroupLayout(jPanel1);
        jPanel1.setLayout(jPanel1Layout);
        jPanel1Layout.setHorizontalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING)
                    .addComponent(emailLabel1, javax.swing.GroupLayout.Alignment.LEADING, 0, 0, Short.MAX_VALUE)
                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                        .addComponent(abinitBannerLabel, javax.swing.GroupLayout.DEFAULT_SIZE, 347, Short.MAX_VALUE)
                        .addComponent(versionLabel)
                        .addComponent(abinitVersionLabel)
                        .addComponent(osLabel)
                        .addComponent(javaVersionLabel)
                        .addComponent(authorLabel)
                        .addComponent(emailLabel)))
                .addContainerGap())
        );
        jPanel1Layout.setVerticalGroup(
            jPanel1Layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(jPanel1Layout.createSequentialGroup()
                .addComponent(abinitBannerLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(versionLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(abinitVersionLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(osLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(javaVersionLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(authorLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.UNRELATED)
                .addComponent(emailLabel)
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                .addComponent(emailLabel1, javax.swing.GroupLayout.DEFAULT_SIZE, 54, Short.MAX_VALUE)
                .addContainerGap())
        );

        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addComponent(jPanel1, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
        );

        pack();
    }// </editor-fold>//GEN-END:initComponents

    // Variables declaration - do not modify//GEN-BEGIN:variables
    private javax.swing.JLabel abinitBannerLabel;
    private javax.swing.JLabel abinitVersionLabel;
    private javax.swing.JLabel authorLabel;
    private javax.swing.JLabel emailLabel;
    private javax.swing.JLabel emailLabel1;
    private javax.swing.JPanel jPanel1;
    private javax.swing.JLabel javaVersionLabel;
    private javax.swing.JLabel osLabel;
    private javax.swing.JLabel versionLabel;
    // End of variables declaration//GEN-END:variables

}
