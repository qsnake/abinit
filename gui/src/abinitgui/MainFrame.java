/*--
MainFrame.java - Created July 16, 2009

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

import java.awt.Color;
import java.io.BufferedInputStream;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.Enumeration;
import java.util.Iterator;
import java.util.List;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.swing.JFileChooser;
import javax.swing.JRadioButton;
import javax.swing.JTable;
import javax.swing.table.DefaultTableColumnModel;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;
import org.jdom.*;

//@SuppressWarnings({"deprecation", "serial"})
public class MainFrame extends javax.swing.JFrame {

    private SFTP sftp = null;
    private SSH ssh = null;
    private MyTableModel pspModel = null;
    private LocalExec localExec;
    private RemoteExec remoteExec;
    private SSHTunnel sshtun;
    private DecimalFormat df_rprim = new DecimalFormat("#0.000000");
    private DecimalFormat df_angdeg = new DecimalFormat("#0.00");
    private AboutDialog about;
    private int lport = 0;
    private DisplayerJDialog outDialog;
    private DisplayerJDialog inputFileDisplayer;
    private DisplayerJDialog clustepInputFileDisplayer;
    private DisplayerJDialog clustepPositionFileDisplayer;
    private GeomDialog geomD;
    private AlCoDialog alcoD;
    private ReReDialog rereD;
    private WaDeDialog wadeD;
    private InOuDialog inouD;
    private TheoDialog theoD;
    //private VarsHelp varsHelp;
    // Avant la version 6 de ABINIT, il y avait deux exécutables différents
    private String SequAbinit = "abinit";
    private String ParaAbinit = "abinit";
    private String CharSet = "UTF-8";

    /** Creates new form MainFrame */
    public MainFrame() {
        initComponents();

        //JFrame.setDefaultLookAndFeelDecorated(true);
        //JDialog.setDefaultLookAndFeelDecorated(true);

        this.setTitle("ABINIT GUI (v. 0.1 April 2011)");

        outDialog = new DisplayerJDialog(this, false);
        outDialog.setTitle("..:: Global MSG Display ::..");

        inputFileDisplayer = new DisplayerJDialog(this, false);
        inputFileDisplayer.setTitle("..:: Input file preview ::..");

        clustepInputFileDisplayer = new DisplayerJDialog(this, false);
        clustepInputFileDisplayer.setTitle("..:: Clustep input file preview ::..");

        clustepPositionFileDisplayer = new DisplayerJDialog(this, false);
        clustepPositionFileDisplayer.setTitle("..:: Clustep position file preview ::..");

        geomD = new GeomDialog(this, false, outDialog);
        geomD.setTitle("..:: Geometry ::..");

        alcoD = new AlCoDialog(this, false, outDialog);
        alcoD.setTitle("..:: Algorithm and convergence ::..");

        rereD = new ReReDialog(this, false, outDialog);
        rereD.setTitle("..:: Real and reciprocal space ::..");

        wadeD = new WaDeDialog(this, false, outDialog);
        wadeD.setTitle("..:: Wavefunctions and densities ::..");

        inouD = new InOuDialog(this, false, outDialog);
        inouD.setTitle("..:: Input / Output ::..");

        theoD = new TheoDialog(this, false, outDialog);
        theoD.setTitle("..:: Theory (DFT) ::..");

        //varsHelp = new VarsHelp();
        //varsHelp.setTitle("..:: Abinit variables help ::..");

        pspModel = new MyTableModel(pspTable);
        pspModel.setNotEditableCol("1-3");
        pspTable.setModel(pspModel);
        initTableHeader(pspTable, new String[]{"Atom", "PSP filename", "PSP type", "PSP path"},
                new Integer[]{null, null, null, null});
        pspTable.setDefaultRenderer(Atom.class,
                new pspAtomRenderer());
        pspTable.setDefaultEditor(Atom.class,
                new AtomEditor(this));

        inputFileTabbedPane.setSelectedIndex(inputFileTabbedPane.getTabCount() - 1);

        localExec = new LocalExec();
        localExec.setDialog(outDialog);
        abinitParaLabel.setText("(local max=" + Runtime.getRuntime().availableProcessors() + ")");

        DecimalFormatSymbols dfs = new DecimalFormatSymbols();
        dfs.setDecimalSeparator('.');
        df_rprim.setDecimalFormatSymbols(dfs);
        df_angdeg.setDecimalFormatSymbols(dfs);

        about = new AboutDialog(this, true);

        loadConfig("config.xml");

        outDialog.setLocationRelativeTo(this);
        outDialog.setVisible(true);

        // TODO rendre visible
        //mainTabbedPane.setEnabledAt(4, false);
    }

    public String getNtypat() {
        return geomD.getNtypat();
    }

    private void initTableHeader(JTable table, String header[], Integer headerWidths[]) {
        TableColumnModel tcm = new DefaultTableColumnModel();
        for (int i = 0; i < header.length; i++) {
            TableColumn tc = new TableColumn(i);
            tc.setHeaderValue(header[i]);
            tc.setResizable(false);
            if (headerWidths[i] != null) {
                tc.setMinWidth(headerWidths[i]);
                tc.setPreferredWidth(headerWidths[i]);
                tc.setMaxWidth(headerWidths[i]);
            }
            tcm.addColumn(tc);
        }
        table.setColumnModel(tcm);
    }

    public void printERR(String s) {
        // TODO mettre de la couleur
        if (s.endsWith("\n")) {
            outDialog.appendERR(s);
        } else {
            outDialog.appendERR(s + "\n");
        }
    }

    public void printOUT(String s) {
        if (s.endsWith("\n")) {
            outDialog.appendOUT(s);
        } else {
            outDialog.appendOUT(s + "\n");
        }
    }

    String removeEndl(String str) {
        if (str.endsWith("\n")) {
            return (String) str.subSequence(0, str.lastIndexOf('\n'));
        } else {
            return str;
        }
    }

    public void printDEB(String str) {
        if (str.endsWith("\n")) {
            outDialog.appendDEB("DEBUG: " + str);
        } else {
            outDialog.appendDEB("DEBUG: " + str + "\n");
        }
    }

    public void printOUT2(String str) {
        if (str.endsWith("\n")) {
            outDialog.appendOUT2(str);
        } else {
            outDialog.appendOUT2(str + "\n");
        }
    }

    public void printGEN(String str, Color color, boolean underline, boolean bolt) {
        if (str.endsWith("\n")) {
            outDialog.appendGEN(str, color, underline, bolt);
        } else {
            outDialog.appendGEN(str + "\n", color, underline, bolt);
        }
    }

    /** This method is called from within the constructor to
     * initialize the form.
     * WARNING: Do NOT modify this code. The content of this method is
     * always regenerated by the Form Editor.
     */
    //@SuppressWarnings("unchecked")
    // <editor-fold defaultstate="collapsed" desc="Generated Code">//GEN-BEGIN:initComponents
    private void initComponents() {

        whereIsAbinitbuttonGroup = new javax.swing.ButtonGroup();
        inputFilebuttonGroup = new javax.swing.ButtonGroup();
        lookAndFeelbuttonGroup = new javax.swing.ButtonGroup();
        abinixbuttonGroup = new javax.swing.ButtonGroup();
        mainTabbedPane = new javax.swing.JTabbedPane();
        configPanel = new javax.swing.JPanel();
        localAbinitRadioButton = new javax.swing.JRadioButton();
        remoteAbinitRadioButton = new javax.swing.JRadioButton();
        whereIsAbinitLabel = new javax.swing.JLabel();
        remoteGatewayRadioButton = new javax.swing.JRadioButton();
        loginPanel = new javax.swing.JPanel();
        hostTextField = new javax.swing.JTextField();
        hostLabel = new javax.swing.JLabel();
        loginLabel = new javax.swing.JLabel();
        loginTextField = new javax.swing.JTextField();
        pwdPasswordField = new javax.swing.JPasswordField();
        pwdLabel = new javax.swing.JLabel();
        gatewayLoginPanel = new javax.swing.JPanel();
        gatewayHostTextField = new javax.swing.JTextField();
        hostBFELabel = new javax.swing.JLabel();
        loginBFELabel = new javax.swing.JLabel();
        gatewayLoginTextField = new javax.swing.JTextField();
        gatewayPasswordField = new javax.swing.JPasswordField();
        pwdBFELabel = new javax.swing.JLabel();
        mySimulationsTextField = new javax.swing.JTextField();
        mySimulationsLabel = new javax.swing.JLabel();
        pspPathTextField = new javax.swing.JTextField();
        pspPathLabel = new javax.swing.JLabel();
        needSGECheckBox = new javax.swing.JCheckBox();
        connectionToggleButton = new javax.swing.JToggleButton();
        abinitPathTextField = new javax.swing.JTextField();
        abinitPathPathLabel = new javax.swing.JLabel();
        sequentialCheckBox = new javax.swing.JCheckBox();
        parallelCheckBox = new javax.swing.JCheckBox();
        abinitParaTextField = new javax.swing.JTextField();
        abinitParaLabel = new javax.swing.JLabel();
        SGEconfigPanel = new javax.swing.JPanel();
        timeLabel = new javax.swing.JLabel();
        nodesLabel = new javax.swing.JLabel();
        ramLabel = new javax.swing.JLabel();
        hdmLabel = new javax.swing.JLabel();
        timeTextField = new javax.swing.JTextField();
        nodesTextField = new javax.swing.JTextField();
        ramTextField = new javax.swing.JTextField();
        hdmTextField = new javax.swing.JTextField();
        emailLabel = new javax.swing.JLabel();
        emailTextField = new javax.swing.JTextField();
        abinitPathButton = new javax.swing.JButton();
        inputFilePanel = new javax.swing.JPanel();
        useCreIFRadioButton = new javax.swing.JRadioButton();
        inputFileTabbedPane = new javax.swing.JTabbedPane();
        basicsScrollPane = new javax.swing.JScrollPane();
        basicsPanel = new javax.swing.JPanel();
        geometryButton = new javax.swing.JButton();
        algoAndConvButton = new javax.swing.JButton();
        realAndRecipButton = new javax.swing.JButton();
        wavefuncAndDensButton = new javax.swing.JButton();
        inputOutputButton = new javax.swing.JButton();
        theoryButton = new javax.swing.JButton();
        jScrollPane5 = new javax.swing.JScrollPane();
        otherTextArea = new javax.swing.JTextArea();
        emptyPanel = new javax.swing.JPanel();
        jLabel4 = new javax.swing.JLabel();
        createButton = new javax.swing.JButton();
        useExtIFRadioButton = new javax.swing.JRadioButton();
        openFileTextField = new javax.swing.JTextField();
        openFileDialogButton = new javax.swing.JButton();
        openFileLabel = new javax.swing.JLabel();
        openXMLFileDialogButton = new javax.swing.JButton();
        openXMLFileTextField = new javax.swing.JTextField();
        openXMLFileLabel = new javax.swing.JLabel();
        saveFileAsButton = new javax.swing.JButton();
        saveFileButton = new javax.swing.JButton();
        sendSIMButton = new javax.swing.JButton();
        pspTextField = new javax.swing.JTextField();
        pspTableScrollPane = new javax.swing.JScrollPane();
        pspTable = new javax.swing.JTable();
        pspLabel = new javax.swing.JLabel();
        displayFileButton = new javax.swing.JButton();
        geditButton = new javax.swing.JButton();
        sshPanel = new javax.swing.JPanel();
        sendSSHButton = new javax.swing.JButton();
        SSHCommandLine = new javax.swing.JTextField();
        startSSHButton = new javax.swing.JButton();
        SSHCommandLineLabel = new javax.swing.JLabel();
        SSHOutputLabel = new javax.swing.JLabel();
        stopSSHButton = new javax.swing.JButton();
        SFTPLoginPanel1 = new javax.swing.JPanel();
        SSHUserAndHostTextField = new javax.swing.JTextField();
        SSHUserAndHostLabel = new javax.swing.JLabel();
        SSHPwdLabel = new javax.swing.JLabel();
        SSHPwdPasswordField = new javax.swing.JPasswordField();
        SSHOutputScrollPane = new javax.swing.JScrollPane();
        SSHOutput = new javax.swing.JTextArea();
        useGlobalConfigSSHCheckBox = new javax.swing.JCheckBox();
        sftpPanel = new javax.swing.JPanel();
        startSFTPButton = new javax.swing.JButton();
        SFTPLoginPanel = new javax.swing.JPanel();
        SFTPUserAndHostTextField = new javax.swing.JTextField();
        SFTPUserAndHostLabel = new javax.swing.JLabel();
        SFTPPwdLabel = new javax.swing.JLabel();
        SFTPPwdPasswordField = new javax.swing.JPasswordField();
        stopSFTPButton = new javax.swing.JButton();
        SFTPCommandLineLabel = new javax.swing.JLabel();
        sendSFTPButton = new javax.swing.JButton();
        SFTPCommandLine = new javax.swing.JTextField();
        SFTPOutputLabel = new javax.swing.JLabel();
        SFTPOutputScrollPane = new javax.swing.JScrollPane();
        SFTPOutput = new javax.swing.JTextArea();
        SFTPProgressBar = new javax.swing.JProgressBar();
        UploadDownloadRATELabel = new javax.swing.JLabel();
        UploadDownloadINFOLabel = new javax.swing.JLabel();
        sendAFileButton = new javax.swing.JButton();
        useGlobalConfigSFTPCheckBox = new javax.swing.JCheckBox();
        mainMenuBar = new javax.swing.JMenuBar();
        fileMenu = new javax.swing.JMenu();
        saveMenuItem = new javax.swing.JMenuItem();
        saveAsMenuItem = new javax.swing.JMenuItem();
        LoadMenuItem = new javax.swing.JMenuItem();
        viewMenu = new javax.swing.JMenu();
        outputMSGMenuItem = new javax.swing.JMenuItem();
        clearOutMSGMenuItem = new javax.swing.JMenuItem();
        postProcMenu = new javax.swing.JMenu();
        getOutputFileMenuItem = new javax.swing.JMenuItem();
        getLogFileMenuItem = new javax.swing.JMenuItem();
        helpMenu = new javax.swing.JMenu();
        helpMenuItem = new javax.swing.JMenuItem();
        varsHelpMenuItem = new javax.swing.JMenuItem();
        jSeparator1 = new javax.swing.JSeparator();
        aboutMenuItem = new javax.swing.JMenuItem();

        setDefaultCloseOperation(javax.swing.WindowConstants.EXIT_ON_CLOSE);
        setTitle("ABINIT GUI (v. 0.1 April 2011)");
        setBackground(new java.awt.Color(245, 242, 239));

        mainTabbedPane.setMaximumSize(new java.awt.Dimension(800, 650));
        mainTabbedPane.setMinimumSize(new java.awt.Dimension(800, 650));
        mainTabbedPane.setPreferredSize(new java.awt.Dimension(800, 650));

        configPanel.setMaximumSize(null);

        whereIsAbinitbuttonGroup.add(localAbinitRadioButton);
        localAbinitRadioButton.setForeground(java.awt.Color.blue);
        localAbinitRadioButton.setText("Local (only for Linux hosts)");
        localAbinitRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                localAbinitRadioButtonActionPerformed(evt);
            }
        });

        whereIsAbinitbuttonGroup.add(remoteAbinitRadioButton);
        remoteAbinitRadioButton.setForeground(java.awt.Color.blue);
        remoteAbinitRadioButton.setText("Remote");
        remoteAbinitRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                remoteAbinitRadioButtonActionPerformed(evt);
            }
        });

        whereIsAbinitLabel.setForeground(java.awt.Color.red);
        whereIsAbinitLabel.setText("ABINIT host location ?");

        whereIsAbinitbuttonGroup.add(remoteGatewayRadioButton);
        remoteGatewayRadioButton.setForeground(java.awt.Color.red);
        remoteGatewayRadioButton.setSelected(true);
        remoteGatewayRadioButton.setText("Remote (behind a gateway)");
        remoteGatewayRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                remoteGatewayRadioButtonActionPerformed(evt);
            }
        });

        loginPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)), "Remote Abinithost login", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("Arial", 3, 14), java.awt.Color.darkGray)); // NOI18N

        hostLabel.setText("hostname or IP");

        loginLabel.setText("Login");

        pwdLabel.setText("Password");

        org.jdesktop.layout.GroupLayout loginPanelLayout = new org.jdesktop.layout.GroupLayout(loginPanel);
        loginPanel.setLayout(loginPanelLayout);
        loginPanelLayout.setHorizontalGroup(
            loginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(loginPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(loginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                    .add(hostLabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(hostTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 125, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(loginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                    .add(loginTextField)
                    .add(loginLabel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 73, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(loginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                    .add(pwdLabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(pwdPasswordField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 100, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        loginPanelLayout.setVerticalGroup(
            loginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(loginPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(loginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(loginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                        .add(hostLabel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 15, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .add(loginLabel))
                    .add(pwdLabel))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(loginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(hostTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(loginTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(pwdPasswordField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        gatewayLoginPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)), "Gateway login", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("Arial", 3, 14), java.awt.Color.darkGray)); // NOI18N

        hostBFELabel.setText("hostname or IP");

        loginBFELabel.setText("Login");

        pwdBFELabel.setText("Password");

        org.jdesktop.layout.GroupLayout gatewayLoginPanelLayout = new org.jdesktop.layout.GroupLayout(gatewayLoginPanel);
        gatewayLoginPanel.setLayout(gatewayLoginPanelLayout);
        gatewayLoginPanelLayout.setHorizontalGroup(
            gatewayLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(gatewayLoginPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(gatewayLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                    .add(hostBFELabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(gatewayHostTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 125, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(gatewayLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                    .add(gatewayLoginTextField)
                    .add(loginBFELabel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 73, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(gatewayLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                    .add(pwdBFELabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(gatewayPasswordField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 100, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        gatewayLoginPanelLayout.setVerticalGroup(
            gatewayLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(gatewayLoginPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(gatewayLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(gatewayLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                        .add(hostBFELabel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 15, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .add(loginBFELabel))
                    .add(pwdBFELabel))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(gatewayLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(gatewayHostTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(gatewayLoginTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(gatewayPasswordField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        mySimulationsLabel.setLabelFor(mySimulationsTextField);
        mySimulationsLabel.setText("Path where to create the simulations filetree");
        mySimulationsLabel.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                mySimulationsLabelMouseClicked(evt);
            }
        });

        pspPathLabel.setLabelFor(pspPathTextField);
        pspPathLabel.setText("Local path to pseudopotentials");
        pspPathLabel.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                pspPathLabelMouseClicked(evt);
            }
        });

        needSGECheckBox.setForeground(new java.awt.Color(0, 128, 0));
        needSGECheckBox.setSelected(true);
        needSGECheckBox.setText("Need SGE script");
        needSGECheckBox.setToolTipText("If checked, a SGE script is created and is used to submit the calculation job.");
        needSGECheckBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                needSGECheckBoxActionPerformed(evt);
            }
        });

        connectionToggleButton.setText("Connect");
        connectionToggleButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                connectionToggleButtonActionPerformed(evt);
            }
        });

        abinitPathPathLabel.setLabelFor(pspPathTextField);
        abinitPathPathLabel.setText("Path to the abinit program (At abinit server !)");
        abinitPathPathLabel.addMouseListener(new java.awt.event.MouseAdapter() {
            public void mouseClicked(java.awt.event.MouseEvent evt) {
                abinitPathPathLabelMouseClicked(evt);
            }
        });

        abinixbuttonGroup.add(sequentialCheckBox);
        sequentialCheckBox.setForeground(java.awt.Color.red);
        sequentialCheckBox.setText("Sequential");
        sequentialCheckBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                sequentialCheckBoxActionPerformed(evt);
            }
        });

        abinixbuttonGroup.add(parallelCheckBox);
        parallelCheckBox.setForeground(java.awt.Color.blue);
        parallelCheckBox.setText("Parallel");
        parallelCheckBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                parallelCheckBoxActionPerformed(evt);
            }
        });

        abinitParaTextField.setEnabled(false);

        abinitParaLabel.setText("(local max = ?)");
        abinitParaLabel.setEnabled(false);

        SGEconfigPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)), "SGE script configuration", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION, javax.swing.border.TitledBorder.DEFAULT_POSITION, new java.awt.Font("DejaVu Sans", 3, 14), java.awt.Color.darkGray)); // NOI18N

        timeLabel.setText("Time [h]");

        nodesLabel.setText("# nodes");

        ramLabel.setText("RAM [Mb]");

        hdmLabel.setText("HDM [Mb]");

        nodesTextField.setEnabled(false);

        hdmTextField.setEnabled(false);

        emailLabel.setText("E-mail where to send feedback");

        org.jdesktop.layout.GroupLayout SGEconfigPanelLayout = new org.jdesktop.layout.GroupLayout(SGEconfigPanel);
        SGEconfigPanel.setLayout(SGEconfigPanelLayout);
        SGEconfigPanelLayout.setHorizontalGroup(
            SGEconfigPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(SGEconfigPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(SGEconfigPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(emailTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 286, Short.MAX_VALUE)
                    .add(SGEconfigPanelLayout.createSequentialGroup()
                        .add(SGEconfigPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING, false)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, timeTextField)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, timeLabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .add(18, 18, 18)
                        .add(SGEconfigPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                            .add(nodesTextField)
                            .add(nodesLabel))
                        .add(18, 18, 18)
                        .add(SGEconfigPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                            .add(ramTextField)
                            .add(ramLabel))
                        .add(18, 18, 18)
                        .add(SGEconfigPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                            .add(hdmTextField)
                            .add(hdmLabel)))
                    .add(emailLabel))
                .addContainerGap())
        );
        SGEconfigPanelLayout.setVerticalGroup(
            SGEconfigPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(SGEconfigPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(SGEconfigPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(SGEconfigPanelLayout.createSequentialGroup()
                        .add(hdmLabel)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(hdmTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                    .add(SGEconfigPanelLayout.createSequentialGroup()
                        .add(ramLabel)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(ramTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                    .add(SGEconfigPanelLayout.createSequentialGroup()
                        .add(nodesLabel)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(nodesTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                    .add(SGEconfigPanelLayout.createSequentialGroup()
                        .add(timeLabel)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(timeTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                .add(emailLabel)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(emailTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );

        abinitPathButton.setText("Find path");
        abinitPathButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                abinitPathButtonActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout configPanelLayout = new org.jdesktop.layout.GroupLayout(configPanel);
        configPanel.setLayout(configPanelLayout);
        configPanelLayout.setHorizontalGroup(
            configPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(configPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(configPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(whereIsAbinitLabel)
                    .add(configPanelLayout.createSequentialGroup()
                        .add(localAbinitRadioButton)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(remoteAbinitRadioButton)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(remoteGatewayRadioButton)
                        .add(18, 18, 18)
                        .add(needSGECheckBox))
                    .add(configPanelLayout.createSequentialGroup()
                        .add(loginPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                        .add(gatewayLoginPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                    .add(mySimulationsLabel)
                    .add(pspPathLabel)
                    .add(abinitPathPathLabel)
                    .add(configPanelLayout.createSequentialGroup()
                        .add(sequentialCheckBox)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(parallelCheckBox)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(abinitParaTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 40, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(abinitParaLabel))
                    .add(SGEconfigPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(configPanelLayout.createSequentialGroup()
                        .add(configPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING, false)
                            .add(abinitPathTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 508, Short.MAX_VALUE)
                            .add(mySimulationsTextField)
                            .add(pspPathTextField))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(configPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(connectionToggleButton, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 143, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                            .add(abinitPathButton))))
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        configPanelLayout.setVerticalGroup(
            configPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(configPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(whereIsAbinitLabel)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                .add(configPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(remoteGatewayRadioButton)
                    .add(localAbinitRadioButton)
                    .add(remoteAbinitRadioButton)
                    .add(needSGECheckBox, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                .add(configPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(configPanelLayout.createSequentialGroup()
                        .add(loginPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .add(18, 18, 18)
                        .add(mySimulationsLabel))
                    .add(gatewayLoginPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(configPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(mySimulationsTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(connectionToggleButton))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(pspPathLabel)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(pspPathTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.UNRELATED)
                .add(abinitPathPathLabel)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(configPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(abinitPathTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(abinitPathButton))
                .add(10, 10, 10)
                .add(configPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(sequentialCheckBox)
                    .add(parallelCheckBox)
                    .add(abinitParaTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(abinitParaLabel))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(SGEconfigPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        mainTabbedPane.addTab("Configuration", configPanel);

        inputFilePanel.setMaximumSize(new java.awt.Dimension(800, 600));
        inputFilePanel.setMinimumSize(new java.awt.Dimension(800, 600));
        inputFilePanel.setPreferredSize(new java.awt.Dimension(800, 600));

        inputFilebuttonGroup.add(useCreIFRadioButton);
        useCreIFRadioButton.setForeground(java.awt.Color.blue);
        useCreIFRadioButton.setText("Use created input file");
        useCreIFRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                useCreIFRadioButtonActionPerformed(evt);
            }
        });

        inputFileTabbedPane.setEnabled(false);
        inputFileTabbedPane.setMaximumSize(new java.awt.Dimension(380, 550));
        inputFileTabbedPane.setMinimumSize(new java.awt.Dimension(380, 550));
        inputFileTabbedPane.setPreferredSize(new java.awt.Dimension(380, 550));

        basicsScrollPane.setEnabled(false);
        basicsScrollPane.setMaximumSize(new java.awt.Dimension(352, 600));
        basicsScrollPane.setMinimumSize(new java.awt.Dimension(352, 600));
        basicsScrollPane.setPreferredSize(new java.awt.Dimension(352, 600));

        basicsPanel.setEnabled(false);
        basicsPanel.setMaximumSize(new java.awt.Dimension(241, 267));
        basicsPanel.setMinimumSize(new java.awt.Dimension(241, 267));

        geometryButton.setText("Geometry");
        geometryButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                geometryButtonActionPerformed(evt);
            }
        });

        algoAndConvButton.setText("Algorithm & convergence");
        algoAndConvButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                algoAndConvButtonActionPerformed(evt);
            }
        });

        realAndRecipButton.setText("Real & reciprocal space");
        realAndRecipButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                realAndRecipButtonActionPerformed(evt);
            }
        });

        wavefuncAndDensButton.setText("Wavefunctions & densities");
        wavefuncAndDensButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                wavefuncAndDensButtonActionPerformed(evt);
            }
        });

        inputOutputButton.setText("Input / Output");
        inputOutputButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                inputOutputButtonActionPerformed(evt);
            }
        });

        theoryButton.setText("Theory (DFT)");
        theoryButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                theoryButtonActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout basicsPanelLayout = new org.jdesktop.layout.GroupLayout(basicsPanel);
        basicsPanel.setLayout(basicsPanelLayout);
        basicsPanelLayout.setHorizontalGroup(
            basicsPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(basicsPanelLayout.createSequentialGroup()
                .add(49, 49, 49)
                .add(basicsPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING, false)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, theoryButton, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, geometryButton, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, algoAndConvButton, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, realAndRecipButton, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, inputOutputButton, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, wavefuncAndDensButton))
                .addContainerGap(91, Short.MAX_VALUE))
        );
        basicsPanelLayout.setVerticalGroup(
            basicsPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(basicsPanelLayout.createSequentialGroup()
                .add(51, 51, 51)
                .add(geometryButton)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(algoAndConvButton)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(realAndRecipButton)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(wavefuncAndDensButton)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(inputOutputButton)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(theoryButton)
                .addContainerGap(125, Short.MAX_VALUE))
        );

        basicsScrollPane.setViewportView(basicsPanel);

        inputFileTabbedPane.addTab("Basics", basicsScrollPane);

        otherTextArea.setColumns(20);
        otherTextArea.setRows(5);
        jScrollPane5.setViewportView(otherTextArea);

        inputFileTabbedPane.addTab("For other variables", jScrollPane5);

        jLabel4.setHorizontalAlignment(javax.swing.SwingConstants.CENTER);
        jLabel4.setText("<HTML> <center> Select the <b>Use created input file</b><br>radiobutton to create graphicaly an input<br>file to send to an ABINIT host.  </HTML>");

        org.jdesktop.layout.GroupLayout emptyPanelLayout = new org.jdesktop.layout.GroupLayout(emptyPanel);
        emptyPanel.setLayout(emptyPanelLayout);
        emptyPanelLayout.setHorizontalGroup(
            emptyPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(emptyPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(jLabel4, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 273, Short.MAX_VALUE)
                .addContainerGap())
        );
        emptyPanelLayout.setVerticalGroup(
            emptyPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(emptyPanelLayout.createSequentialGroup()
                .add(133, 133, 133)
                .add(jLabel4, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap(203, Short.MAX_VALUE))
        );

        inputFileTabbedPane.addTab("", emptyPanel);

        createButton.setText("Generate file preview (test)");
        createButton.setEnabled(false);
        createButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                createButtonActionPerformed(evt);
            }
        });

        inputFilebuttonGroup.add(useExtIFRadioButton);
        useExtIFRadioButton.setForeground(java.awt.Color.red);
        useExtIFRadioButton.setSelected(true);
        useExtIFRadioButton.setText("Use an external input file");
        useExtIFRadioButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                useExtIFRadioButtonActionPerformed(evt);
            }
        });

        openFileDialogButton.setText("...");
        openFileDialogButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                openFileDialogButtonActionPerformed(evt);
            }
        });

        openFileLabel.setText("Open the ABINIT input file (usualy *.in)");

        openXMLFileDialogButton.setText("...");
        openXMLFileDialogButton.setEnabled(false);
        openXMLFileDialogButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                openXMLFileDialogButtonActionPerformed(evt);
            }
        });

        openXMLFileTextField.setEnabled(false);

        openXMLFileLabel.setText("Open the ABINIT input file (usualy *.ab [XML file])");
        openXMLFileLabel.setEnabled(false);

        saveFileAsButton.setText("Save As");
        saveFileAsButton.setEnabled(false);

        saveFileButton.setText("Save");
        saveFileButton.setEnabled(false);

        sendSIMButton.setText("<HTML> <center> <b>Send the simulation</b><br> the simulation will start at server side </HTML>");
        sendSIMButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                sendSIMButtonActionPerformed(evt);
            }
        });

        pspTextField.addKeyListener(new java.awt.event.KeyAdapter() {
            public void keyReleased(java.awt.event.KeyEvent evt) {
                pspTextFieldKeyReleased(evt);
            }
        });

        pspTable.setModel(new javax.swing.table.DefaultTableModel(
            new Object [][] {

            },
            new String [] {

            }
        ));
        pspTable.setRowSelectionAllowed(false);
        pspTable.getTableHeader().setReorderingAllowed(false);
        pspTableScrollPane.setViewportView(pspTable);

        pspLabel.setText("Number of different psps (=ntypat)");

        displayFileButton.setText("Display");
        displayFileButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                displayFileButtonActionPerformed(evt);
            }
        });

        geditButton.setText("Edit");
        geditButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                geditButtonActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout inputFilePanelLayout = new org.jdesktop.layout.GroupLayout(inputFilePanel);
        inputFilePanel.setLayout(inputFilePanelLayout);
        inputFilePanelLayout.setHorizontalGroup(
            inputFilePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(inputFilePanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(inputFilePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(inputFilePanelLayout.createSequentialGroup()
                        .add(pspLabel)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(pspTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 54, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                    .add(org.jdesktop.layout.GroupLayout.TRAILING, inputFilePanelLayout.createSequentialGroup()
                        .add(inputFilePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, sendSIMButton, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 424, Short.MAX_VALUE)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, useExtIFRadioButton)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, inputFilePanelLayout.createSequentialGroup()
                                .add(openFileLabel)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(displayFileButton)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(geditButton))
                            .add(org.jdesktop.layout.GroupLayout.LEADING, inputFilePanelLayout.createSequentialGroup()
                                .add(openFileDialogButton)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(openFileTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 340, Short.MAX_VALUE))
                            .add(org.jdesktop.layout.GroupLayout.LEADING, openXMLFileLabel)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, useCreIFRadioButton)
                            .add(org.jdesktop.layout.GroupLayout.LEADING, inputFilePanelLayout.createSequentialGroup()
                                .add(saveFileButton)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(saveFileAsButton)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(createButton))
                            .add(org.jdesktop.layout.GroupLayout.LEADING, inputFilePanelLayout.createSequentialGroup()
                                .add(openXMLFileDialogButton)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(openXMLFileTextField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 340, Short.MAX_VALUE))
                            .add(org.jdesktop.layout.GroupLayout.LEADING, pspTableScrollPane, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 424, Short.MAX_VALUE))
                        .add(14, 14, 14)))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(inputFileTabbedPane, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 334, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );
        inputFilePanelLayout.setVerticalGroup(
            inputFilePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(inputFilePanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(inputFilePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING, false)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, inputFileTabbedPane, 0, 0, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, inputFilePanelLayout.createSequentialGroup()
                        .add(useExtIFRadioButton)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(inputFilePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(openFileLabel)
                            .add(displayFileButton)
                            .add(geditButton))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(inputFilePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(openFileDialogButton)
                            .add(openFileTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(inputFilePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(pspLabel)
                            .add(pspTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(pspTableScrollPane, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 98, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .add(17, 17, 17)
                        .add(useCreIFRadioButton)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(openXMLFileLabel)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(inputFilePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(openXMLFileDialogButton)
                            .add(openXMLFileTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(inputFilePanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(saveFileButton)
                            .add(saveFileAsButton)
                            .add(createButton))
                        .add(18, 18, 18)
                        .add(sendSIMButton, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)))
                .add(145, 145, 145))
        );

        inputFilePanelLayout.linkSize(new java.awt.Component[] {openFileDialogButton, openFileTextField}, org.jdesktop.layout.GroupLayout.VERTICAL);

        inputFilePanelLayout.linkSize(new java.awt.Component[] {openXMLFileDialogButton, openXMLFileTextField}, org.jdesktop.layout.GroupLayout.VERTICAL);

        mainTabbedPane.addTab("Input File", inputFilePanel);

        sshPanel.setBackground(new java.awt.Color(186, 195, 199));
        sshPanel.setMaximumSize(new java.awt.Dimension(800, 600));
        sshPanel.setMinimumSize(new java.awt.Dimension(800, 600));
        sshPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        sendSSHButton.setText("Send");
        sendSSHButton.setEnabled(false);
        sendSSHButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                sendSSHButtonActionPerformed(evt);
            }
        });

        SSHCommandLine.setEnabled(false);
        SSHCommandLine.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                SSHCommandLineActionPerformed(evt);
            }
        });

        startSSHButton.setText("Start SSH client");
        startSSHButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                startSSHButtonActionPerformed(evt);
            }
        });

        SSHCommandLineLabel.setText("Command line");

        SSHOutputLabel.setText("Output");

        stopSSHButton.setText("Stop SSH client");
        stopSSHButton.setEnabled(false);
        stopSSHButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                stopSSHButtonActionPerformed(evt);
            }
        });

        SFTPLoginPanel1.setBackground(new java.awt.Color(186, 195, 199));

        SSHUserAndHostTextField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                SSHUserAndHostTextFieldActionPerformed(evt);
            }
        });

        SSHUserAndHostLabel.setText("Username@Hostname");

        SSHPwdLabel.setText("Password");

        SSHPwdPasswordField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                SSHPwdPasswordFieldActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout SFTPLoginPanel1Layout = new org.jdesktop.layout.GroupLayout(SFTPLoginPanel1);
        SFTPLoginPanel1.setLayout(SFTPLoginPanel1Layout);
        SFTPLoginPanel1Layout.setHorizontalGroup(
            SFTPLoginPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(SFTPLoginPanel1Layout.createSequentialGroup()
                .addContainerGap()
                .add(SFTPLoginPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(SSHUserAndHostLabel)
                    .add(SSHUserAndHostTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 196, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(SFTPLoginPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(SSHPwdLabel)
                    .add(SSHPwdPasswordField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 83, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
        );
        SFTPLoginPanel1Layout.setVerticalGroup(
            SFTPLoginPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(SFTPLoginPanel1Layout.createSequentialGroup()
                .add(SFTPLoginPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(SSHUserAndHostLabel)
                    .add(SSHPwdLabel))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(SFTPLoginPanel1Layout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(SSHUserAndHostTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(SSHPwdPasswordField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );

        SSHOutputScrollPane.setAutoscrolls(true);

        SSHOutput.setColumns(20);
        SSHOutput.setEditable(false);
        SSHOutput.setFont(new java.awt.Font("Courier New", 0, 12));
        SSHOutput.setRows(5);
        SSHOutput.setEnabled(false);
        SSHOutputScrollPane.setViewportView(SSHOutput);

        useGlobalConfigSSHCheckBox.setBackground(new java.awt.Color(186, 195, 199));
        useGlobalConfigSSHCheckBox.setText("Use the global configuration");
        useGlobalConfigSSHCheckBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                useGlobalConfigSSHCheckBoxActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout sshPanelLayout = new org.jdesktop.layout.GroupLayout(sshPanel);
        sshPanel.setLayout(sshPanelLayout);
        sshPanelLayout.setHorizontalGroup(
            sshPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(org.jdesktop.layout.GroupLayout.TRAILING, sshPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(sshPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.TRAILING)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, SSHOutputScrollPane, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 693, Short.MAX_VALUE)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, sshPanelLayout.createSequentialGroup()
                        .add(sshPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(sshPanelLayout.createSequentialGroup()
                                .add(startSSHButton)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(stopSSHButton))
                            .add(useGlobalConfigSSHCheckBox))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(SFTPLoginPanel1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                    .add(org.jdesktop.layout.GroupLayout.LEADING, SSHCommandLineLabel)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, SSHOutputLabel)
                    .add(org.jdesktop.layout.GroupLayout.LEADING, sshPanelLayout.createSequentialGroup()
                        .add(SSHCommandLine, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 609, Short.MAX_VALUE)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(sendSSHButton)))
                .addContainerGap())
        );
        sshPanelLayout.setVerticalGroup(
            sshPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(sshPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(sshPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(sshPanelLayout.createSequentialGroup()
                        .add(sshPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                            .add(startSSHButton)
                            .add(stopSSHButton))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(useGlobalConfigSSHCheckBox)
                        .add(3, 3, 3)
                        .add(SSHCommandLineLabel))
                    .add(SFTPLoginPanel1, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(sshPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(SSHCommandLine, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(sendSSHButton))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(SSHOutputLabel)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(SSHOutputScrollPane, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 380, Short.MAX_VALUE)
                .addContainerGap())
        );

        mainTabbedPane.addTab("SSH terminal", sshPanel);

        sftpPanel.setBackground(new java.awt.Color(186, 195, 199));
        sftpPanel.setMaximumSize(new java.awt.Dimension(800, 600));
        sftpPanel.setMinimumSize(new java.awt.Dimension(800, 600));
        sftpPanel.setPreferredSize(new java.awt.Dimension(800, 600));

        startSFTPButton.setText("Start SFTP client");
        startSFTPButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                startSFTPButtonActionPerformed(evt);
            }
        });

        SFTPLoginPanel.setBackground(new java.awt.Color(186, 195, 199));
        SFTPLoginPanel.setPreferredSize(new java.awt.Dimension(321, 58));

        SFTPUserAndHostTextField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                SFTPUserAndHostTextFieldActionPerformed(evt);
            }
        });

        SFTPUserAndHostLabel.setText("Username@Hostname");

        SFTPPwdLabel.setText("Password");

        SFTPPwdPasswordField.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                SFTPPwdPasswordFieldActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout SFTPLoginPanelLayout = new org.jdesktop.layout.GroupLayout(SFTPLoginPanel);
        SFTPLoginPanel.setLayout(SFTPLoginPanelLayout);
        SFTPLoginPanelLayout.setHorizontalGroup(
            SFTPLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(SFTPLoginPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(SFTPLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(SFTPUserAndHostLabel)
                    .add(SFTPUserAndHostTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 191, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(SFTPLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(SFTPPwdLabel)
                    .add(SFTPPwdPasswordField, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 80, Short.MAX_VALUE))
                .addContainerGap())
        );
        SFTPLoginPanelLayout.setVerticalGroup(
            SFTPLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(SFTPLoginPanelLayout.createSequentialGroup()
                .add(SFTPLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(SFTPUserAndHostLabel)
                    .add(SFTPPwdLabel))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(SFTPLoginPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(SFTPUserAndHostTextField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(SFTPPwdPasswordField, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );

        stopSFTPButton.setText("Stop SFTP client");
        stopSFTPButton.setEnabled(false);
        stopSFTPButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                stopSFTPButtonActionPerformed(evt);
            }
        });

        SFTPCommandLineLabel.setText("Command line");

        sendSFTPButton.setText("Send");
        sendSFTPButton.setEnabled(false);
        sendSFTPButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                sendSFTPButtonActionPerformed(evt);
            }
        });

        SFTPCommandLine.setEnabled(false);
        SFTPCommandLine.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                SFTPCommandLineActionPerformed(evt);
            }
        });

        SFTPOutputLabel.setText("Output");

        SFTPOutputScrollPane.setAutoscrolls(true);

        SFTPOutput.setColumns(20);
        SFTPOutput.setEditable(false);
        SFTPOutput.setFont(new java.awt.Font("Courier New", 0, 12));
        SFTPOutput.setRows(5);
        SFTPOutput.setEnabled(false);
        SFTPOutputScrollPane.setViewportView(SFTPOutput);

        UploadDownloadRATELabel.setText("Upload / Download RATE");

        UploadDownloadINFOLabel.setText("Upload / Download INFO");

        sendAFileButton.setText("Send a file");
        sendAFileButton.setEnabled(false);
        sendAFileButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                sendAFileButtonActionPerformed(evt);
            }
        });

        useGlobalConfigSFTPCheckBox.setBackground(new java.awt.Color(186, 195, 199));
        useGlobalConfigSFTPCheckBox.setText("Use the global configuration");
        useGlobalConfigSFTPCheckBox.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                useGlobalConfigSFTPCheckBoxActionPerformed(evt);
            }
        });

        org.jdesktop.layout.GroupLayout sftpPanelLayout = new org.jdesktop.layout.GroupLayout(sftpPanel);
        sftpPanel.setLayout(sftpPanelLayout);
        sftpPanelLayout.setHorizontalGroup(
            sftpPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(sftpPanelLayout.createSequentialGroup()
                .addContainerGap()
                .add(sftpPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(SFTPOutputScrollPane, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 752, Short.MAX_VALUE)
                    .add(SFTPOutputLabel)
                    .add(org.jdesktop.layout.GroupLayout.TRAILING, sftpPanelLayout.createSequentialGroup()
                        .add(SFTPCommandLine, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 668, Short.MAX_VALUE)
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(sendSFTPButton))
                    .add(SFTPCommandLineLabel)
                    .add(UploadDownloadINFOLabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 752, Short.MAX_VALUE)
                    .add(UploadDownloadRATELabel, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 752, Short.MAX_VALUE)
                    .add(SFTPProgressBar, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 752, Short.MAX_VALUE)
                    .add(sftpPanelLayout.createSequentialGroup()
                        .add(sftpPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(sftpPanelLayout.createSequentialGroup()
                                .add(startSFTPButton)
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(stopSFTPButton))
                            .add(useGlobalConfigSFTPCheckBox))
                        .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                        .add(SFTPLoginPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                        .add(18, 18, 18)
                        .add(sendAFileButton)))
                .addContainerGap())
        );
        sftpPanelLayout.setVerticalGroup(
            sftpPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(sftpPanelLayout.createSequentialGroup()
                .add(sftpPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                    .add(sftpPanelLayout.createSequentialGroup()
                        .addContainerGap()
                        .add(sftpPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
                            .add(sftpPanelLayout.createSequentialGroup()
                                .add(sftpPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                                    .add(startSFTPButton)
                                    .add(stopSFTPButton))
                                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                                .add(useGlobalConfigSFTPCheckBox)
                                .add(3, 3, 3)
                                .add(SFTPCommandLineLabel))
                            .add(SFTPLoginPanel, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 58, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)))
                    .add(sftpPanelLayout.createSequentialGroup()
                        .add(31, 31, 31)
                        .add(sendAFileButton)))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(sftpPanelLayout.createParallelGroup(org.jdesktop.layout.GroupLayout.BASELINE)
                    .add(SFTPCommandLine, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                    .add(sendSFTPButton))
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(SFTPOutputLabel)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(SFTPOutputScrollPane, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, 304, Short.MAX_VALUE)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(UploadDownloadINFOLabel)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(UploadDownloadRATELabel)
                .addPreferredGap(org.jdesktop.layout.LayoutStyle.RELATED)
                .add(SFTPProgressBar, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        mainTabbedPane.addTab("SFTP terminal", sftpPanel);

        fileMenu.setLabel("File");

        saveMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_S, java.awt.event.InputEvent.CTRL_MASK));
        saveMenuItem.setText("Save config.");
        saveMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                saveMenuItemActionPerformed(evt);
            }
        });
        fileMenu.add(saveMenuItem);

        saveAsMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_S, java.awt.event.InputEvent.SHIFT_MASK | java.awt.event.InputEvent.CTRL_MASK));
        saveAsMenuItem.setText("Save config. as...");
        saveAsMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                saveAsMenuItemActionPerformed(evt);
            }
        });
        fileMenu.add(saveAsMenuItem);

        LoadMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_L, java.awt.event.InputEvent.CTRL_MASK));
        LoadMenuItem.setText("Load config.");
        LoadMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                LoadMenuItemActionPerformed(evt);
            }
        });
        fileMenu.add(LoadMenuItem);

        mainMenuBar.add(fileMenu);

        viewMenu.setText("View");

        outputMSGMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_D, java.awt.event.InputEvent.CTRL_MASK));
        outputMSGMenuItem.setText("Display MSG panel");
        outputMSGMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                outputMSGMenuItemActionPerformed(evt);
            }
        });
        viewMenu.add(outputMSGMenuItem);

        clearOutMSGMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_R, java.awt.event.InputEvent.CTRL_MASK));
        clearOutMSGMenuItem.setText("Clear MSG panel");
        clearOutMSGMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                clearOutMSGMenuItemActionPerformed(evt);
            }
        });
        viewMenu.add(clearOutMSGMenuItem);

        mainMenuBar.add(viewMenu);

        postProcMenu.setText("PostProc");

        getOutputFileMenuItem.setText("Download output file");
        getOutputFileMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                getOutputFileMenuItemActionPerformed(evt);
            }
        });
        postProcMenu.add(getOutputFileMenuItem);

        getLogFileMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_L, java.awt.event.InputEvent.ALT_MASK));
        getLogFileMenuItem.setText("Download log file");
        getLogFileMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                getLogFileMenuItemActionPerformed(evt);
            }
        });
        postProcMenu.add(getLogFileMenuItem);

        mainMenuBar.add(postProcMenu);

        helpMenu.setMnemonic('h');
        helpMenu.setText("Help");

        helpMenuItem.setText("Help contents");
        helpMenuItem.setEnabled(false);
        helpMenu.add(helpMenuItem);

        varsHelpMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_H, java.awt.event.InputEvent.ALT_MASK));
        varsHelpMenuItem.setText("Abinit variables help");
        varsHelpMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                varsHelpMenuItemActionPerformed(evt);
            }
        });
        helpMenu.add(varsHelpMenuItem);
        helpMenu.add(jSeparator1);

        aboutMenuItem.setAccelerator(javax.swing.KeyStroke.getKeyStroke(java.awt.event.KeyEvent.VK_A, java.awt.event.InputEvent.ALT_MASK));
        aboutMenuItem.setText("About");
        aboutMenuItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                aboutMenuItemActionPerformed(evt);
            }
        });
        helpMenu.add(aboutMenuItem);

        mainMenuBar.add(helpMenu);

        setJMenuBar(mainMenuBar);

        org.jdesktop.layout.GroupLayout layout = new org.jdesktop.layout.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(org.jdesktop.layout.GroupLayout.TRAILING, layout.createSequentialGroup()
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .add(mainTabbedPane, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 754, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(org.jdesktop.layout.GroupLayout.LEADING)
            .add(org.jdesktop.layout.GroupLayout.TRAILING, layout.createSequentialGroup()
                .addContainerGap(org.jdesktop.layout.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .add(mainTabbedPane, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE, 613, org.jdesktop.layout.GroupLayout.PREFERRED_SIZE)
                .addContainerGap())
        );

        getAccessibleContext().setAccessibleName("");

        pack();
    }// </editor-fold>//GEN-END:initComponents

    private void stopConnection() {
        if (remoteGatewayRadioButton.isSelected()) {
            if (remoteExec != null) {
                remoteExec.stop();
                remoteExec = null;
            }
            if (sshtun != null) {
                sshtun.stop();
                sshtun = null;
            }
            lport = 0;
        } else if (remoteAbinitRadioButton.isSelected()) {
            if (remoteExec != null) {
                remoteExec.stop();
                remoteExec = null;
            }
        } else if (localAbinitRadioButton.isSelected()) {
            // Pas besoin en local
        } else { // Le choix n'a pas été fait
            printERR("Choose a destination option please at config. tab !");
        }
    }

    private void aboutMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_aboutMenuItemActionPerformed
        about.setLocationRelativeTo(this);
        about.setVisible(true);
    }//GEN-LAST:event_aboutMenuItemActionPerformed

    private void saveMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_saveMenuItemActionPerformed
        saveConfig(null);
    }//GEN-LAST:event_saveMenuItemActionPerformed

    private void LoadMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_LoadMenuItemActionPerformed
        JFileChooser fc = new JFileChooser(".");
        fc.setMultiSelectionEnabled(false);
        int retValue = fc.showOpenDialog(this);
        if (retValue == JFileChooser.APPROVE_OPTION) {
            File file = fc.getSelectedFile();
            loadConfig(file.getPath());
        }

    }//GEN-LAST:event_LoadMenuItemActionPerformed

    private void saveAsMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_saveAsMenuItemActionPerformed
        JFileChooser fc = new JFileChooser(".");
        fc.setMultiSelectionEnabled(false);
        int retValue = fc.showSaveDialog(this);
        if (retValue == JFileChooser.APPROVE_OPTION) {
            File file = fc.getSelectedFile();
            saveConfig(file.getPath());
        }
    }//GEN-LAST:event_saveAsMenuItemActionPerformed

    private void outputMSGMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_outputMSGMenuItemActionPerformed
        outDialog.setLocationRelativeTo(this);
        outDialog.setVisible(true);
    }//GEN-LAST:event_outputMSGMenuItemActionPerformed

    private void clearOutMSGMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_clearOutMSGMenuItemActionPerformed
        outDialog.clear();
    }//GEN-LAST:event_clearOutMSGMenuItemActionPerformed

    private void getOutputFileMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_getOutputFileMenuItemActionPerformed
        Runnable r = new Runnable() {

            @Override
            public void run() {
                if ((remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) && remoteExec == null) {
                    printERR("Please connect to a ABINIT host before downloading anything!");
                    sendSIMButton.setEnabled(true);
                    return;
                }

                String rootPath = mySimulationsTextField.getText();
                String outputFolder = "output";

                String inputFile = "";
                String inputFileName = "";

                if (useExtIFRadioButton.isSelected()) {
                    inputFile = openFileTextField.getText();
                    inputFileName = Utils.getLastToken(inputFile.replace('\\', '/'), "/");
                } else if (useCreIFRadioButton.isSelected()) {
                    inputFile = openXMLFileTextField.getText();
                    inputFileName = Utils.getLastToken(inputFile.replace('\\', '/'), "/");
                } else {
                    printERR("Choose an option please ! (use an external inputfile or created a inputfile)");
                    sendSIMButton.setEnabled(true);
                    return;
                }

                // Test de l'existance de inputfile
                if (!Utils.exists(inputFile)) {
                    printERR("The file " + inputFile + " doesn't exist !");
                    sendSIMButton.setEnabled(true);
                    return;
                }

                String simName = null;
                if (inputFileName != null) {
                    if (!inputFileName.equals("")) {
                        int idx = inputFileName.indexOf('.');
                        if (idx > 0 && idx < inputFileName.length()) {
                            simName = inputFileName.substring(0, idx);
                        } else {
                            simName = inputFileName;
                        }
                    } else {
                        printERR("inputFileName == \"\"");
                        return;
                    }
                } else {
                    printERR("inputFileName == null");
                    return;
                }

                if (!inputFile.equals("")) {

                    String outputPath = rootPath + "/" + outputFolder;
                    String fileName = outputPath + "/" + simName + ".out";

                    // Réception (copie) du fichier d'output si celui-ci est distant
                    if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
                        String file = "";
                        String outputFiles = getOutputFilesR(fileName + "*");
                        StringTokenizer st = new StringTokenizer(outputFiles, "\n");
                        while (st.hasMoreElements()) {
                            file = st.nextToken();
                            printOUT2("File = " + file);
                            if (Utils.osName().startsWith("Windows")) {
                                sendCommand("unix2dos " + file);
                            }
                            getFile(file + " " + file);
                            if (Utils.osName().startsWith("Windows")) {
                                sendCommand("dos2unix " + file);
                                // TODO Util.unix2dos(new File(file))
                            }
                        }
                        fileName = file; // Prend le nom du dernier fichier!
                    }

                    // ****************************************************************************
                    // Tester l'existance du fichier
                    if (!Utils.exists(fileName)) {
                        printERR("File " + fileName + " doesn't exist !");
                        return;
                    } else {
                        if (Utils.osName().equals("Linux")) {
                            localCommand("gedit " + fileName);
                        } else if (Utils.osName().equals("Mac OS X")) {
                            localCommand("open -a textedit " + fileName);
                        } else {
                            localCommand("notepad " + fileName);
                        }
                    }
                    // ****************************************************************************
                } else {
                    printERR("Please setup the inputfile textfield !");
                    return;
                }
            }
        };

        Thread t = new Thread(r);
        t.start();
    }//GEN-LAST:event_getOutputFileMenuItemActionPerformed

    private void varsHelpMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_varsHelpMenuItemActionPerformed
        //varsHelp.setVisible(true);
        printERR("Not implemented yet!");
    }//GEN-LAST:event_varsHelpMenuItemActionPerformed

    private void getLogFileMenuItemActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_getLogFileMenuItemActionPerformed
        Runnable r = new Runnable() {

            @Override
            public void run() {
                if ((remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) && remoteExec == null) {
                    printERR("Please connect to a ABINIT host before downloading anything!");
                    sendSIMButton.setEnabled(true);
                    return;
                }

                String rootPath = mySimulationsTextField.getText();
                //String inputFolder = "input";
                //String outputFolder = "output";
                //String wholedataFolder = "wholedata";
                //String pseudopotFolder = "pseudopot";
                String logfilesFolder = "logfiles";

                String inputFile = "";
                String inputFileName = "";

                if (useExtIFRadioButton.isSelected()) {
                    inputFile = openFileTextField.getText();
                    inputFileName = Utils.getLastToken(inputFile.replace('\\', '/'), "/");
                } else if (useCreIFRadioButton.isSelected()) {
                    inputFile = openXMLFileTextField.getText();
                    inputFileName = Utils.getLastToken(inputFile.replace('\\', '/'), "/");
                } else {
                    printERR("Choose an option please ! (use an external inputfile or created a inputfile)");
                    sendSIMButton.setEnabled(true);
                    return;
                }

                // Test de l'existance de inputfile
                if (!Utils.exists(inputFile)) {
                    printERR("The file " + inputFile + " doesn't exist !");
                    sendSIMButton.setEnabled(true);
                    return;
                }

                String simName = null;
                if (inputFileName != null) {
                    if (!inputFileName.equals("")) {
                        int idx = inputFileName.indexOf('.');
                        if (idx > 0 && idx < inputFileName.length()) {
                            simName = inputFileName.substring(0, idx);
                        } else {
                            simName = inputFileName;
                        }
                    } else {
                        printERR("inputFileName == \"\"");
                        return;
                    }
                } else {
                    printERR("inputFileName == null");
                    return;
                }

                if (!inputFile.equals("")) {
                    String fileName = rootPath + "/" + logfilesFolder + "/" + simName + ".log";

                    // Réception (copie) du fichier *.log si celui-ci est distant
                    if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
                        if (Utils.osName().startsWith("Windows")) {
                            sendCommand("unix2dos " + fileName);
                        }
                        getFile(fileName + " " + fileName);
                        if (Utils.osName().startsWith("Windows")) {
                            sendCommand("dos2unix " + fileName);
                            // TODO Util.unix2dos(new File(fileName))
                        }
                    }

                    // ****************************************************************************
                    // Tester l'existance du fichier
                    if (!Utils.exists(fileName)) {
                        printERR("File " + fileName + " doesn't exist !");
                        return;
                    } else {
                        if (Utils.osName().equals("Linux")) {
                            localCommand("gedit " + fileName);
                        } else if (Utils.osName().equals("Mac OS X")) {
                            localCommand("open -a textedit " + fileName);
                        } else {
                            localCommand("notepad " + fileName);
                        }
                    }
                    // ****************************************************************************
                } else {
                    printERR("Please setup the inputfile textfield !");
                    return;
                }
            }
        };

        Thread t = new Thread(r);
        t.start();
    }//GEN-LAST:event_getLogFileMenuItemActionPerformed

    private void useGlobalConfigSFTPCheckBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_useGlobalConfigSFTPCheckBoxActionPerformed
        // TODO add your handling code here:
        if (useGlobalConfigSFTPCheckBox.isSelected()) {
            if (remoteGatewayRadioButton.isSelected()) {
                if (lport != 0) {
                    startSFTPButton.setEnabled(false);
                    SFTPUserAndHostTextField.setEnabled(false);
                    SFTPPwdPasswordField.setEnabled(false);
                    SFTPUserAndHostLabel.setEnabled(false);
                    SFTPPwdLabel.setEnabled(false);

                    if (sftp != null) {
                        sftp.stop();
                    }
                    sftp = new SFTP(SFTPProgressBar, UploadDownloadINFOLabel, UploadDownloadRATELabel, SFTPOutput);
                    sftp.setDialog(outDialog);
                    sftp.setUserAndHost(loginTextField.getText() + "@localhost");
                    sftp.setPassword(new String(pwdPasswordField.getPassword()));
                    sftp.setPort(lport);

                    if (sftp.start()) {
                        sendSFTPButton.setEnabled(true);
                        stopSFTPButton.setEnabled(false);
                        SFTPCommandLine.setEnabled(true);
                        SFTPOutput.setEnabled(true);
                        sendAFileButton.setEnabled(true);
                    } else {
                        useGlobalConfigSFTPCheckBox.setSelected(false);
                        startSFTPButton.setEnabled(true);
                        SFTPUserAndHostTextField.setEnabled(true);
                        SFTPPwdPasswordField.setEnabled(true);
                        SFTPUserAndHostLabel.setEnabled(true);
                        SFTPPwdLabel.setEnabled(true);
                        printERR("The SFTP client is not well connected. Please connect properly"
                                + " to an abinit Host at the configuration tab before !");
                    }
                } else {
                    useGlobalConfigSFTPCheckBox.setSelected(false);
                    printERR("The ssh tunnel is not working, please connect at the config. tab befor !");
                }
            } else if (remoteAbinitRadioButton.isSelected()) {
                startSFTPButton.setEnabled(false);
                SFTPUserAndHostTextField.setEnabled(false);
                SFTPPwdPasswordField.setEnabled(false);
                SFTPUserAndHostLabel.setEnabled(false);
                SFTPPwdLabel.setEnabled(false);

                if (sftp != null) {
                    sftp.stop();
                }
                sftp = new SFTP(SFTPProgressBar, UploadDownloadINFOLabel, UploadDownloadRATELabel, SFTPOutput);
                sftp.setDialog(outDialog);
                sftp.setUserAndHost(loginTextField.getText() + "@" + hostTextField.getText());
                sftp.setPassword(new String(pwdPasswordField.getPassword()));

                if (sftp.start()) {
                    sendSFTPButton.setEnabled(true);
                    stopSFTPButton.setEnabled(false);
                    SFTPCommandLine.setEnabled(true);
                    SFTPOutput.setEnabled(true);
                    sendAFileButton.setEnabled(true);
                } else {
                    useGlobalConfigSFTPCheckBox.setSelected(false);
                    startSFTPButton.setEnabled(true);
                    SFTPUserAndHostTextField.setEnabled(true);
                    SFTPPwdPasswordField.setEnabled(true);
                    SFTPUserAndHostLabel.setEnabled(true);
                    SFTPPwdLabel.setEnabled(true);
                    printERR("The SFTP client is not well connected. Please connect properly"
                            + " to an abinit Host at the configuration tab before !");
                }
            } else if (localAbinitRadioButton.isSelected()) {
                useGlobalConfigSFTPCheckBox.setSelected(false);
                printERR("To connect to the local host, please choose the remote option"
                        + " where you specify localhost as hostname !");
            } else { // Le choix n'a pas été fait
                printERR("Choose a destination option please at config. tab !");
            }
        } else {
            if (sftp != null) {
                sftp.stop();
            }
            startSFTPButton.setEnabled(true);
            SFTPUserAndHostTextField.setEnabled(true);
            SFTPPwdPasswordField.setEnabled(true);
            SFTPUserAndHostLabel.setEnabled(true);
            SFTPPwdLabel.setEnabled(true);
            sendSFTPButton.setEnabled(false);
            SFTPCommandLine.setEnabled(false);
            sendAFileButton.setEnabled(false);
        }
}//GEN-LAST:event_useGlobalConfigSFTPCheckBoxActionPerformed

    private void sendAFileButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_sendAFileButtonActionPerformed
        sendSFTPButton.setEnabled(false);
        SFTPCommandLine.setEditable(false);
        new SendSFTPThread(mySimulationsTextField.getText());
}//GEN-LAST:event_sendAFileButtonActionPerformed

    private void SFTPCommandLineActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_SFTPCommandLineActionPerformed
        sendSFTPButtonActionPerformed(evt);
}//GEN-LAST:event_SFTPCommandLineActionPerformed

    private void sendSFTPButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_sendSFTPButtonActionPerformed
        sendSFTPButton.setEnabled(false);
        SFTPCommandLine.setEditable(false);
        new SFTPCMDThread(SFTPCommandLine.getText());
}//GEN-LAST:event_sendSFTPButtonActionPerformed

    private void stopSFTPButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_stopSFTPButtonActionPerformed
        stopSFTPButton.setEnabled(false);

        if (sftp != null) {
            sftp.stop();
            sftp = null;
        }

        startSFTPButton.setEnabled(true);
        SFTPUserAndHostTextField.setEnabled(true);
        SFTPPwdPasswordField.setEnabled(true);
        SFTPUserAndHostLabel.setEnabled(true);
        SFTPPwdLabel.setEnabled(true);
        SFTPCommandLine.setEnabled(false);
        SFTPOutput.setEnabled(false);
        sendSFTPButton.setEnabled(false);
        sendAFileButton.setEnabled(false);
}//GEN-LAST:event_stopSFTPButtonActionPerformed

    private void SFTPPwdPasswordFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_SFTPPwdPasswordFieldActionPerformed
        startSFTPButtonActionPerformed(evt);
}//GEN-LAST:event_SFTPPwdPasswordFieldActionPerformed

    private void SFTPUserAndHostTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_SFTPUserAndHostTextFieldActionPerformed
        startSFTPButtonActionPerformed(evt);
}//GEN-LAST:event_SFTPUserAndHostTextFieldActionPerformed

    private void startSFTPButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_startSFTPButtonActionPerformed
        startSFTPButton.setEnabled(false);
        SFTPUserAndHostTextField.setEnabled(false);
        SFTPPwdPasswordField.setEnabled(false);
        SFTPUserAndHostLabel.setEnabled(false);
        SFTPPwdLabel.setEnabled(false);
        //startSFTPButton.setBackground(new Color(0, 255, 0));

        sftp = new SFTP(SFTPProgressBar, UploadDownloadINFOLabel, UploadDownloadRATELabel, SFTPOutput);

        sftp.setUserAndHost(SFTPUserAndHostTextField.getText());
        sftp.setPassword(new String(SFTPPwdPasswordField.getPassword()));

        if (sftp.start()) {
            sendSFTPButton.setEnabled(true);
            stopSFTPButton.setEnabled(true);
            SFTPCommandLine.setEnabled(true);
            SFTPOutput.setEnabled(true);
            sendAFileButton.setEnabled(true);
        } else {
            startSFTPButton.setEnabled(true);
            SFTPUserAndHostTextField.setEnabled(true);
            SFTPPwdPasswordField.setEnabled(true);
            SFTPUserAndHostLabel.setEnabled(true);
            SFTPPwdLabel.setEnabled(true);
        }
}//GEN-LAST:event_startSFTPButtonActionPerformed

    private void useGlobalConfigSSHCheckBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_useGlobalConfigSSHCheckBoxActionPerformed
        // TODO add your handling code here:
        if (useGlobalConfigSSHCheckBox.isSelected()) {
            if (remoteGatewayRadioButton.isSelected()) {
                if (lport != 0) {
                    startSSHButton.setEnabled(false);
                    SSHUserAndHostTextField.setEnabled(false);
                    SSHPwdPasswordField.setEnabled(false);
                    SSHUserAndHostLabel.setEnabled(false);
                    SSHPwdLabel.setEnabled(false);

                    if (ssh != null) {
                        ssh.stop();
                    }
                    ssh = new SSH(SSHOutput, true);

                    ssh.setUserAndHost(loginTextField.getText() + "@localhost");
                    ssh.setPassword(new String(pwdPasswordField.getPassword()));
                    ssh.setPort(lport);

                    if (ssh.start()) {
                        sendSSHButton.setEnabled(true);
                        stopSSHButton.setEnabled(false);
                        SSHCommandLine.setEnabled(true);
                        SSHOutput.setEnabled(true);
                    } else {
                        useGlobalConfigSSHCheckBox.setSelected(false);
                        startSSHButton.setEnabled(true);
                        SSHUserAndHostTextField.setEnabled(true);
                        SSHPwdPasswordField.setEnabled(true);
                        SSHUserAndHostLabel.setEnabled(true);
                        SSHPwdLabel.setEnabled(true);
                        printERR("The SSH client is not well connected. Please connect properly"
                                + " to an abinit Host at the configuration tab before !");
                    }
                } else {
                    useGlobalConfigSSHCheckBox.setSelected(false);
                    printERR("The ssh tunnel is not working, please connect at the config. tab befor !");
                }
            } else if (remoteAbinitRadioButton.isSelected()) {
                startSSHButton.setEnabled(false);
                SSHUserAndHostTextField.setEnabled(false);
                SSHPwdPasswordField.setEnabled(false);
                SSHUserAndHostLabel.setEnabled(false);
                SSHPwdLabel.setEnabled(false);

                if (ssh != null) {
                    ssh.stop();
                }
                ssh = new SSH(SSHOutput, true);
                ssh.setDialog(outDialog);
                ssh.setUserAndHost(loginTextField.getText() + "@" + hostTextField.getText());
                ssh.setPassword(new String(pwdPasswordField.getPassword()));

                if (ssh.start()) {
                    sendSSHButton.setEnabled(true);
                    stopSSHButton.setEnabled(false);
                    SSHCommandLine.setEnabled(true);
                    SSHOutput.setEnabled(true);
                } else {
                    useGlobalConfigSSHCheckBox.setSelected(false);
                    startSSHButton.setEnabled(true);
                    SSHUserAndHostTextField.setEnabled(true);
                    SSHPwdPasswordField.setEnabled(true);
                    SSHUserAndHostLabel.setEnabled(true);
                    SSHPwdLabel.setEnabled(true);
                    printERR("The SSH client is not well connected. Please connect properly"
                            + " to an abinit Host at the configuration tab before !");
                }
            } else if (localAbinitRadioButton.isSelected()) {
                useGlobalConfigSSHCheckBox.setSelected(false);
                printERR("To connect to the local host, please choose the remote option"
                        + " where you specify localhost as hostname !");
            } else { // Le choix n'a pas été fait
                printERR("Choose a destination option please at config. tab !");
            }
        } else {
            if (ssh != null) {
                ssh.stop();
            }
            startSSHButton.setEnabled(true);
            SSHUserAndHostTextField.setEnabled(true);
            SSHPwdPasswordField.setEnabled(true);
            SSHUserAndHostLabel.setEnabled(true);
            SSHPwdLabel.setEnabled(true);
            sendSSHButton.setEnabled(false);
            SSHCommandLine.setEnabled(false);
            SSHOutput.setText("");
        }
}//GEN-LAST:event_useGlobalConfigSSHCheckBoxActionPerformed

    private void SSHPwdPasswordFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_SSHPwdPasswordFieldActionPerformed
        startSSHButtonActionPerformed(evt);
}//GEN-LAST:event_SSHPwdPasswordFieldActionPerformed

    private void SSHUserAndHostTextFieldActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_SSHUserAndHostTextFieldActionPerformed
        startSSHButtonActionPerformed(evt);
}//GEN-LAST:event_SSHUserAndHostTextFieldActionPerformed

    private void stopSSHButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_stopSSHButtonActionPerformed
        stopSSHButton.setEnabled(false);

        if (ssh != null) {
            ssh.stop();
        }

        startSSHButton.setEnabled(true);
        SSHUserAndHostTextField.setEnabled(true);
        SSHPwdPasswordField.setEnabled(true);
        SSHUserAndHostLabel.setEnabled(true);
        SSHPwdLabel.setEnabled(true);
        SSHCommandLine.setEnabled(false);
        SSHOutput.setEnabled(false);
        sendSSHButton.setEnabled(false);
}//GEN-LAST:event_stopSSHButtonActionPerformed

    private void startSSHButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_startSSHButtonActionPerformed
        startSSHButton.setEnabled(false);
        SSHUserAndHostTextField.setEnabled(false);
        SSHPwdPasswordField.setEnabled(false);
        SSHUserAndHostLabel.setEnabled(false);
        SSHPwdLabel.setEnabled(false);

        ssh = new SSH(SSHOutput, true);

        ssh.setUserAndHost(SSHUserAndHostTextField.getText());
        ssh.setPassword(new String(SSHPwdPasswordField.getPassword()));

        if (ssh.start()) {
            sendSSHButton.setEnabled(true);
            stopSSHButton.setEnabled(true);
            SSHCommandLine.setEnabled(true);
            SSHOutput.setEnabled(true);
        } else {
            startSSHButton.setEnabled(true);
            SSHUserAndHostTextField.setEnabled(true);
            SSHPwdPasswordField.setEnabled(true);
            SSHUserAndHostLabel.setEnabled(true);
            SSHPwdLabel.setEnabled(true);
        }
}//GEN-LAST:event_startSSHButtonActionPerformed

    private void SSHCommandLineActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_SSHCommandLineActionPerformed
        sendSSHButtonActionPerformed(evt);
}//GEN-LAST:event_SSHCommandLineActionPerformed

    private void sendSSHButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_sendSSHButtonActionPerformed
        ssh.sendCommand(SSHCommandLine.getText());
        SSHCommandLine.setText("");
}//GEN-LAST:event_sendSSHButtonActionPerformed

    private void geditButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_geditButtonActionPerformed
        Runnable r = new Runnable() {

            @Override
            public void run() {
                String fileName = openFileTextField.getText();

                // ****************************************************************************
                // Tester l'existance du fichier
                if (!Utils.exists(fileName)) {
                    printERR("File " + fileName + " doesn't exist !");
                    return;
                } else {
                    if (Utils.osName().equals("Linux")) {
                        localCommand("gedit " + fileName);
                    } else if (Utils.osName().equals("Mac OS X")) {
                        localCommand("open -a textedit " + fileName);
                    } else if (Utils.osName().startsWith("Windows")) {
                        localCommand("notepad " + fileName);
                    } else {
                        printDEB("You must be in a UNIX platform to edit the input"
                                + "\nfile with a text editor from this GUI!");
                    }
                }
                // ****************************************************************************
            }
        };

        Thread t = new Thread(r);
        t.start();
}//GEN-LAST:event_geditButtonActionPerformed

    private void displayFileButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_displayFileButtonActionPerformed
        inputFileDisplayer.setVisible(true);
        // TODO : pour quand ce sera éditable
        //inputFileDisplayer.setEditable(true);

        String fileContent = "";

        try {
            File file = new File(openFileTextField.getText());

            /*FileReader fr = new FileReader(file);
            BufferedReader br = new BufferedReader(fr);

            String thisLine;
            while ((thisLine = br.readLine()) != null) {
                fileContent += thisLine + "\n";
            }*/

            Scanner scanner = new Scanner(file).useDelimiter("\\Z");
            fileContent = scanner.next();
            scanner.close();

        } catch (FileNotFoundException e) {
            printERR(e.getMessage());
        } catch (IOException e) {
            printERR(e.getMessage());
        }

        inputFileDisplayer.setText(fileContent);
}//GEN-LAST:event_displayFileButtonActionPerformed

    private void pspTextFieldKeyReleased(java.awt.event.KeyEvent evt) {//GEN-FIRST:event_pspTextFieldKeyReleased
        try {
            int npsp = Integer.parseInt(pspTextField.getText());

            if (npsp > 1000) {
                npsp = 1000;
                Object strTab[][] = new Object[npsp][4];
                for (int i = 0; i < npsp; i++) {
                    strTab[i] = new Object[]{new Atom(), "", "", ""};
                }
                pspModel.setData(strTab);
                //znuclTable.setModel(znuclModel);
            } else {
                Object strTab[][] = new Object[npsp][4];
                for (int i = 0; i < npsp; i++) {
                    strTab[i] = new Object[]{new Atom(), "", "", ""};
                }
                pspModel.setData(strTab);
                //znuclTable.setModel(znuclModel);
            }
        } catch (Exception e) {
            //printERR(e.getMessage());
            pspModel.setData(null);
            //znuclTable.setModel(znuclModel);
        }
}//GEN-LAST:event_pspTextFieldKeyReleased

    private void sendSIMButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_sendSIMButtonActionPerformed
        sendSIMButton.setEnabled(false);
        Runnable r = new Runnable() {

            @Override
            public void run() {
                if (localAbinitRadioButton.isSelected() && Utils.osName().startsWith("Windows")) {
                    printERR("Please connect to a remote CLUSTEP host before submitting a simulation !");
                    sendSIMButton.setEnabled(true);
                    return;
                }

                if ((remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) && remoteExec == null) {
                    printERR("Please connect to a ABINIT host before submitting a simulation !");
                    sendSIMButton.setEnabled(true);
                    return;
                }

                createFiletree();

                String rootPath = mySimulationsTextField.getText();
                String inputFolder = "input";
                String outputFolder = "output";
                String wholedataFolder = "wholedata";
                String pseudopotFolder = "pseudopot";
                String logfilesFolder = "logfiles";

                // ********************************************************************************************************************************

                String cwd = "";

                String CMD = "pwd";

                RetMSG retmsg;
                if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
                    if (remoteExec != null) {
                        retmsg = remoteExec.sendCommand(CMD);
                        if (retmsg.getRetCode() == RetMSG.SUCCES) {
                            printOUT("PWD: " + retmsg.getRetMSG());
                            cwd = removeEndl(retmsg.getRetMSG());
                        } else {
                            //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
                            printERR("Error: " + retmsg.getRetMSG() + " !");
                        }
                    } else {
                        printERR("First connect to an abinit host please !");
                    }
                } else if (localAbinitRadioButton.isSelected()) {
                    if (localExec != null) {
                        retmsg = localExec.sendCommand(CMD);
                        if (retmsg.getRetCode() == RetMSG.SUCCES) {
                            printOUT("PWD: " + retmsg.getRetMSG());
                            cwd = removeEndl(retmsg.getRetMSG());
                        } else {
                            //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
                            printERR("Error: " + retmsg.getRetMSG() + " !");
                        }
                    }
                } else { // Le choix n'a pas été fait
                    printERR("Choose a destination option please at config. tab !");
                }

                // ********************************************************************************************************************************

                String inputFile = "";
                String inputFileName = "";
                String pathToAbinit = abinitPathTextField.getText();

                String sep = Utils.fileSeparator();

                if (useExtIFRadioButton.isSelected()) {
                    inputFile = openFileTextField.getText();
                    inputFileName = Utils.getLastToken(inputFile.replace('\\', '/'), "/");
                } else if (useCreIFRadioButton.isSelected()) {
                    inputFile = openXMLFileTextField.getText();
                    // TODO créer le fichier input (*.in)

                    // Création du fichier input
                    try {
                        OutputStreamWriter fw = new OutputStreamWriter(new FileOutputStream(inputFile), CharSet);
                        //FileWriter fw = new FileWriter(inputFile);
                        BufferedWriter bw = new BufferedWriter(fw);
                        PrintWriter pw = new PrintWriter(bw);
                        // écriture des variable basics dans le fichier input
                        pw.print(geomD.getData()); // TODO
                        pw.print(alcoD.getData()); // TODO
                        pw.print(rereD.getData()); // TODO
                        pw.print(wadeD.getData()); // TODO
                        pw.print(inouD.getData()); // TODO
                        pw.print(theoD.getData()); // TODO
                        pw.print(otherTextArea.getText()); // TODO
                        pw.close();
                        bw.close();
                        fw.close();
                    } catch (IOException e) {
                        //printERR("Exception in sendSIMButtonActionPerformed:" + e + "");
                        printERR("The input file could not be created !");
                        sendSIMButton.setEnabled(true);
                        return;
                    }

                    inputFileName = Utils.getLastToken(inputFile.replace('\\', '/'), "/");
                } else {
                    printERR("Choose an option please ! (use an external inputfile or created a inputfile)");
                    sendSIMButton.setEnabled(true);
                    return;
                }

                // Test de l'existance de inputfile
                if (!Utils.exists(inputFile)) {
                    printERR("The file " + inputFile + " doesn't exist !");
                    sendSIMButton.setEnabled(true);
                    return;
                }

                String simName = null;
                if (inputFileName != null) {
                    int idx = inputFileName.indexOf('.');
                    if (idx > 0 && idx < inputFileName.length()) {
                        simName = inputFileName.substring(0, idx);
                    } else {
                        simName = inputFileName;
                    }
                }

                if (!inputFile.equals("")) {
                    int nbProc;
                    if (abinitParaTextField.isEnabled()) {
                        try {
                            nbProc = Integer.parseInt(abinitParaTextField.getText());
                        } catch (Exception e) {
                            printERR("Please set up the number of processors to use !");
                            sendSIMButton.setEnabled(true);
                            return;
                        }
                    } else {
                        nbProc = 1;
                    }

                    if (needSGECheckBox.isSelected()) {
                        int time, nodes, ram, hdm;
                        String email;
                        try {
                            time = Integer.parseInt(timeTextField.getText());
                            //nodes = Integer.parseInt(nodesTextField.getText());
                            ram = Integer.parseInt(ramTextField.getText());
                            //hdm = Integer.parseInt(hdmTextField.getText());
                            email = emailTextField.getText();
                        } catch (Exception e) {
                            //printERR("Exception in sendSIMButtonActionPerformed:" + e + "");
                            printERR("The SGE script configurations are probably wrong !");
                            sendSIMButton.setEnabled(true);
                            return;
                        }
                        // Création du script SGE
                        try {
                            //String rootPath_ = (new File(rootPath)).getCanonicalPath();

                            String PBSfileName = rootPath + sep + simName + ".SGE.sh";
                            OutputStreamWriter fw = new OutputStreamWriter(new FileOutputStream(PBSfileName), CharSet);
                            //FileWriter fw = new FileWriter(PBSfileName);
                            BufferedWriter bw = new BufferedWriter(fw);
                            PrintWriter pw = new PrintWriter(bw);

                            //*********************************************************************************************
                            String fileContent = "#!/bin/bash" + "\n"
                                    + "#" + "\n"
                                    + "# On old Green node" + "\n"
                                    + "#$ -l nb=false" + "\n"
                                    + "#" + "\n"
                                    + "# Ask for pe=parrallel environment, snode or openmpi" + "\n"
                                    + "# snode= same node, as the shared memory communication is the fastest" + "\n"
                                    + "#$ -pe openmpi " + nbProc + "\n"
                                    + "# -pe snode8 8" + "\n"
                                    + "\n"
                                    + "# keep current working directory" + "\n"
                                    + "#$ -cwd" + "\n"
                                    + "\n"
                                    + "#$ -o SGE_out-$JOB_ID.log" + "\n"
                                    + "#$ -e SGE_err-$JOB_ID.log" + "\n"
                                    + "\n"
                                    + "# give a name to your job" + "\n"
                                    + "#$ -N " + simName + "\n"
                                    + "\n"
                                    + "# keep all the defined variables" + "\n"
                                    + "#$ -V" + "\n"
                                    + "#$ -l nb=false" + "\n"
                                    + "\n"
                                    + "# not mandatory: highmem=true (hm=true) for 32GB node" + "\n"
                                    + "# or hm=false for 16GB node" + "\n"
                                    + "# no hm argument does not take about the kind of node ram (16/32)" + "\n"
                                    + "# -l hm=true" + "\n"
                                    + "\n"
                                    + "# IMPORTANT: You need to specify the mem_free" + "\n"
                                    + "# h_vmem can also be set but mf is mandatory!" + "\n"
                                    + "# max 31G if hm=true and max 15G if hm=false" + "\n"
                                    + "#$ -l mf=" + ram + "M" + "\n"
                                    + "\n"
                                    + "# Specify the requested time" + "\n"
                                    + "#$ -l h_rt=" + time + ":00:00" + "\n"
                                    + "\n"
                                    + "# To be informed by email (besa= begin,end,stop,abort)" + "\n"
                                    + "#$ -M " + email + "\n"
                                    + "#$ -m besa" + "\n"
                                    //+ "# ---------------------------" + "\n"
                                    + "\n";
                            if (parallelCheckBox.isSelected()) {
                                fileContent += "MPI=/cvos/shared/apps/openmpi/intel/64/1.3.1/bin/mpirun" + "\n"
                                        + "${MPI} -np " + nbProc + " " + pathToAbinit + "/" + ParaAbinit + " < " + cwd + "/"
                                        + rootPath.replaceFirst("./", "") + "/" + simName + ".files >& " + cwd + "/"
                                        + rootPath.replaceFirst("./", "") + "/" + logfilesFolder + "/" + simName + ".log";
                            } else {
                                fileContent += pathToAbinit + "/" + ParaAbinit + " < " + cwd + "/"
                                        + rootPath.replaceFirst("./", "") + "/" + simName + ".files >& " + cwd + "/"
                                        + rootPath.replaceFirst("./", "") + "/" + logfilesFolder + "/" + simName + ".log";
                            }
                            pw.print(fileContent);
                            //*********************************************************************************************

                            pw.println();
                            pw.close();
                            bw.close();
                            fw.close();
                        } catch (IOException e) {
                            //printERR("Exception in sendSIMButtonActionPerformed:" + e + "");
                            printERR("The SGE script could not be created !");
                            sendSIMButton.setEnabled(true);
                            return;
                        }
                    } else {
                        // Création du script BASH
                        try {
                            String SHfileName = rootPath + sep + simName + ".sh";
                            OutputStreamWriter fw = new OutputStreamWriter(new FileOutputStream(SHfileName), CharSet);
                            //FileWriter fw = new FileWriter(SHfileName);
                            BufferedWriter bw = new BufferedWriter(fw);
                            PrintWriter pw = new PrintWriter(bw);
                            pw.println("#!/bin/bash");

                            if (parallelCheckBox.isSelected()) {
                                pw.print("mpirun -np " + nbProc + " " + pathToAbinit + "/" + ParaAbinit + " < " + cwd + "/" + rootPath.replaceFirst("./", "") + "/" + simName + ".files >& "
                                        + cwd + "/" + rootPath.replaceFirst("./", "") + "/" + logfilesFolder + "/" + simName + ".log");
                            } else {
                                pw.print(pathToAbinit + "/" + SequAbinit + " < " + cwd + "/" + rootPath.replaceFirst("./", "") + "/" + simName + ".files >& "
                                        + cwd + "/" + rootPath.replaceFirst("./", "") + "/" + logfilesFolder + "/" + simName + ".log");
                            }

                            pw.println();
                            pw.close();
                            bw.close();
                            fw.close();
                        } catch (IOException e) {
                            //printDEB("Exception in sendSIMButtonActionPerformed:" + e + "");
                            printERR("The bash script could not be created !");
                            sendSIMButton.setEnabled(true);
                            return;
                        }
                    }

                    // Envoie (copie) du fichier d'input
                    String inputFileR = rootPath + "/" + inputFolder + "/" + inputFileName;
                    putFile(inputFile + " " + inputFileR);

                    if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
                        if (Utils.osName().startsWith("Windows")) {
                            sendCommand("dos2unix " + inputFileR);
                            // TODO Util.dos2unix(new File(inputFileR)); // Transformer avant d'envoyer le fichier
                        }
                    }

                    // Création du contenu du fichier de configuration (*.files)
                    String configFileContent = "";
                    configFileContent += cwd + "/" + rootPath.replaceFirst("./", "") + "/" + inputFolder + "/" + inputFileName + "\n";
                    configFileContent += cwd + "/" + rootPath.replaceFirst("./", "") + "/" + outputFolder + "/" + simName + ".out\n";
                    configFileContent += cwd + "/" + rootPath.replaceFirst("./", "") + "/" + wholedataFolder + "/" + simName + "/" + simName + "i\n";
                    configFileContent += cwd + "/" + rootPath.replaceFirst("./", "") + "/" + wholedataFolder + "/" + simName + "/" + simName + "o\n";
                    configFileContent += cwd + "/" + rootPath.replaceFirst("./", "") + "/" + wholedataFolder + "/" + simName + "/" + simName + "\n";

                    if (useExtIFRadioButton.isSelected()) {
                        int row = pspTable.getRowCount();
                        if (row > 0) {
                            for (int i = 0; i < row; i++) {
                                try {
                                    Atom at = (Atom) pspTable.getValueAt(i, 0);
                                    // Envoie du fichier avec le pseudopotentiel
                                    putFile(at.getPSPPath() + sep + at.getPSPFileName()
                                            + " " + rootPath + "/" + pseudopotFolder + "/" + at.getPSPFileName());
                                    configFileContent += cwd + "/" + rootPath.replaceFirst("./", "") + "/" + pseudopotFolder + "/" + at.getPSPFileName() + "\n";
                                } catch (Exception e) {
                                    //printERR(e.getMessage());
                                    printERR("Error processing pseudopotential(s) !");
                                    sendSIMButton.setEnabled(true);
                                    return;
                                }
                            }
                        } else {
                            printERR("Please setup pseudopotential(s) !");
                            sendSIMButton.setEnabled(true);
                            return;
                        }
                    } else if (useCreIFRadioButton.isSelected()) {
                        int row = geomD.getZnuclTable().getRowCount();
                        if (row > 0) {
                            for (int i = 0; i < row; i++) {
                                try {
                                    Atom at = (Atom) geomD.getZnuclTable().getValueAt(i, 0);
                                    // Envoie du fichier avec le pseudopotentiel
                                    putFile(at.getPSPPath() + sep + at.getPSPFileName()
                                            + " " + rootPath + "/" + pseudopotFolder + "/" + at.getPSPFileName());
                                    configFileContent += cwd + "/" + rootPath.replaceFirst("./", "") + "/" + pseudopotFolder + "/" + at.getPSPFileName() + "\n";
                                } catch (Exception e) {
                                    //printERR(e.getMessage());
                                    printERR("Error processing pseudopotential(s) !");
                                    sendSIMButton.setEnabled(true);
                                    return;
                                }
                            }
                        } else {
                            printERR("Please setup pseudopotential(s) !");
                            sendSIMButton.setEnabled(true);
                            return;
                        }
                    } else {
                        printERR("Choose an option please ! (use an external inputfile or created a inputfile)");
                        sendSIMButton.setEnabled(true);
                        return;
                    }

                    // Création du fichier de configuration
                    try {
                        String FILESfileName = rootPath + sep + simName + ".files";
                        OutputStreamWriter fw = new OutputStreamWriter(new FileOutputStream(FILESfileName), CharSet);
                        //FileWriter fw = new FileWriter(FILESfileName);
                        BufferedWriter bw = new BufferedWriter(fw);
                        PrintWriter pw = new PrintWriter(bw);
                        pw.print(configFileContent);
                        pw.close();
                        bw.close();
                        fw.close();
                    } catch (IOException e) {
                        //printERR(e.getMessage());
                        printERR("The configuration file (*.files) could not be created !");
                        sendSIMButton.setEnabled(true);
                        return;
                    }
                    if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
                        // Envoie du fichier de configuration
                        String configFile = rootPath + sep + simName + ".files";
                        String configFileR = rootPath + "/" + simName + ".files";

                        /*if (Utils.osName().startsWith("Windows")) {
                                Utils.dos2unix(new File(configFile));
                            }*/

                        putFile(configFile + " " + configFileR);

                        if (Utils.osName().startsWith("Windows")) {
                            sendCommand("dos2unix " + configFileR);
                        }
                    }

                    // Creation du dossier simName dans wholedataFolder
                    mkdir(rootPath + "/" + wholedataFolder + "/" + simName);
                    if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
                        mkdirR(rootPath + "/" + wholedataFolder + "/" + simName);
                    }

                    if (needSGECheckBox.isSelected()) {
                        String sgeSHFile = rootPath + sep + simName + ".SGE.sh";
                        String sgeSHFileR = rootPath + "/" + simName + ".SGE.sh";
                        if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {

                            /*if (Utils.osName().startsWith("Windows")) {
                                Utils.dos2unix(new File(sgeSHFileR));
                            }*/

                            // Envoie du fichier SGE
                            putFile(sgeSHFile + " " + sgeSHFileR);

                            if (Utils.osName().startsWith("Windows")) {
                                sendCommand("dos2unix " + sgeSHFileR);
                            }
                        }
                        // lancement des commandes d'exécution de la simulation
                        sendCommand("qsub " + sgeSHFileR);
                    } else {
                        String SHFile = rootPath + sep + simName + ".sh";
                        String SHFileR = rootPath + "/" + simName + ".sh";
                        if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
                            /*if (Utils.osName().startsWith("Windows")) {
                                Utils.dos2unix(new File(SHFileR));
                            }*/
                            // Envoie du fichier BASH
                            putFile(SHFile + " " + SHFileR);

                            if (Utils.osName().startsWith("Windows")) {
                                sendCommand("dos2unix " + SHFileR);
                            }
                        }
                        // lancement des commandes d'exécution de la simulation
                        //printDEB("bash " + SHFileR);
                        sendCommand("bash " + SHFileR);
                        //printDEB("bash " + SHFileR);
                    }
                } else {
                    printERR("Please setup the inputfile textfield !");
                    sendSIMButton.setEnabled(true);
                    return;
                }

                if (localAbinitRadioButton.isSelected()) {
                    printOUT("The simulation was submitted to the local Abinit server.");
                } else {
                    printOUT("The simulation was submitted to the remote Abinit server " + hostTextField.getText());
                    if (remoteGatewayRadioButton.isSelected()) {
                        printOUT(" via the gateway " + gatewayHostTextField.getText() + ".");
                    } else {
                        //printOUT(".");
                    }
                }
                printDEB("The submission thread ended successfully! (Abinit)");
                sendSIMButton.setEnabled(true);
            }
        };

        Thread t = new Thread(r);
        t.start();
}//GEN-LAST:event_sendSIMButtonActionPerformed

    private void openXMLFileDialogButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_openXMLFileDialogButtonActionPerformed
        JFileChooser fc = new JFileChooser(mySimulationsTextField.getText());
        fc.setMultiSelectionEnabled(false);
        int retValue = fc.showOpenDialog(this);
        if (retValue == JFileChooser.APPROVE_OPTION) {
            File file = fc.getSelectedFile();
            openXMLFileTextField.setText(file.getPath());
        }
}//GEN-LAST:event_openXMLFileDialogButtonActionPerformed

    private void openFileDialogButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_openFileDialogButtonActionPerformed
        JFileChooser fc = new JFileChooser(".");
        File currDir = new File(".");
        String currPath = currDir.getAbsolutePath();
        String basePath = basePath = currPath.replace("\\", "/").replace(".", "");
        printDEB(basePath);
        fc.setMultiSelectionEnabled(false);

        int retValue = fc.showOpenDialog(this);
        if (retValue == JFileChooser.APPROVE_OPTION) {
            File file = fc.getSelectedFile();
            String relPath = file.getAbsolutePath().replace("\\", "/").replace(basePath, "./");
            openFileTextField.setText(relPath);
        }
}//GEN-LAST:event_openFileDialogButtonActionPerformed

    private void useExtIFRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_useExtIFRadioButtonActionPerformed
        useExtIFRadioButton.setForeground(Color.red);
        useCreIFRadioButton.setForeground(Color.blue);

        openFileLabel.setEnabled(true);
        openFileDialogButton.setEnabled(true);
        displayFileButton.setEnabled(true);
        geditButton.setEnabled(true);
        openFileTextField.setEnabled(true);
        pspLabel.setEnabled(true);
        pspTextField.setEnabled(true);
        pspTable.setEnabled(true);
        pspTable.setVisible(true);

        openXMLFileLabel.setEnabled(false);
        openXMLFileDialogButton.setEnabled(false);
        openXMLFileTextField.setEnabled(false);
        saveFileAsButton.setEnabled(false);
        saveFileButton.setEnabled(false);
        createButton.setEnabled(false);
        inputFileDisplayer.setVisible(false);
        inputFileTabbedPane.setEnabled(false);

        inputFileTabbedPane.setSelectedIndex(inputFileTabbedPane.getTabCount() - 1);
}//GEN-LAST:event_useExtIFRadioButtonActionPerformed

    private void createButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_createButtonActionPerformed
        inputFileDisplayer.setVisible(true);
        //inputFileDisplayer.setText(getBasics()); // TODO
        inputFileDisplayer.setText(geomD.getData() + alcoD.getData() + rereD.getData()
                + wadeD.getData() + inouD.getData() + theoD.getData() + otherTextArea.getText());
}//GEN-LAST:event_createButtonActionPerformed

    private void theoryButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_theoryButtonActionPerformed
        // TODO add your handling code here:
        theoD.setVisible(true);
}//GEN-LAST:event_theoryButtonActionPerformed

    private void inputOutputButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_inputOutputButtonActionPerformed
        // TODO add your handling code here:
        inouD.setVisible(true);
}//GEN-LAST:event_inputOutputButtonActionPerformed

    private void wavefuncAndDensButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_wavefuncAndDensButtonActionPerformed
        // TODO add your handling code here:
        wadeD.setVisible(true);
}//GEN-LAST:event_wavefuncAndDensButtonActionPerformed

    private void realAndRecipButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_realAndRecipButtonActionPerformed
        // TODO add your handling code here:
        rereD.setVisible(true);
}//GEN-LAST:event_realAndRecipButtonActionPerformed

    private void algoAndConvButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_algoAndConvButtonActionPerformed
        // TODO add your handling code here:
        alcoD.setVisible(true);
}//GEN-LAST:event_algoAndConvButtonActionPerformed

    private void geometryButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_geometryButtonActionPerformed
        // TODO add your handling code here:
        geomD.setVisible(true);
}//GEN-LAST:event_geometryButtonActionPerformed

    private void useCreIFRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_useCreIFRadioButtonActionPerformed
        useCreIFRadioButton.setForeground(Color.red);
        useExtIFRadioButton.setForeground(Color.blue);

        openXMLFileLabel.setEnabled(true);
        openXMLFileDialogButton.setEnabled(true);
        openXMLFileTextField.setEnabled(true);
        saveFileAsButton.setEnabled(true);
        saveFileButton.setEnabled(true);
        createButton.setEnabled(true);
        inputFileTabbedPane.setEnabled(true);

        openFileLabel.setEnabled(false);
        openFileDialogButton.setEnabled(false);
        displayFileButton.setEnabled(false);
        geditButton.setEnabled(false);
        openFileTextField.setEnabled(false);
        pspLabel.setEnabled(false);
        pspTextField.setEnabled(false);
        pspTable.setEnabled(false);
        pspTable.setVisible(false);

        inputFileTabbedPane.setSelectedIndex(0);
}//GEN-LAST:event_useCreIFRadioButtonActionPerformed

    private void abinitPathButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_abinitPathButtonActionPerformed
        abinitPathButton.setEnabled(false);
        Runnable r = new Runnable() {

            @Override
            public void run() {
                String CMD = "whereis " + SequAbinit;
                RetMSG retmsg;
                if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
                    if (remoteExec != null) {
                        retmsg = remoteExec.sendCommand(CMD);
                        if (retmsg.getRetCode() == RetMSG.SUCCES) {
                            StringTokenizer st = new StringTokenizer(retmsg.getRetMSG());
                            int nbt = st.countTokens();
                            for (int i = 0; i < nbt; i++) {
                                String str = st.nextToken();
                                if (i == 1) {
                                    // TODO adapter au systÃ¨me Windows
                                    int idx = str.lastIndexOf('/');
                                    abinitPathTextField.setText((String) str.subSequence(0, idx));
                                }
                                printOUT(str);
                            }
                        } else {
                            //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
                            printERR("Error: " + retmsg.getRetMSG() + " !");
                        }
                    } else {
                        printERR("First connect to an abinit host please !");
                    }
                } else if (localAbinitRadioButton.isSelected()) {
                    if (localExec != null) {
                        retmsg = localExec.sendCommand(CMD);
                        if (retmsg.getRetCode() == RetMSG.SUCCES) {
                            StringTokenizer st = new StringTokenizer(retmsg.getRetMSG());
                            int nbt = st.countTokens();
                            for (int i = 0; i < nbt; i++) {
                                String str = st.nextToken();
                                if (i == 1) {
                                    // TODO adapter au systÃ¨me Windows
                                    int idx = str.lastIndexOf('/');
                                    abinitPathTextField.setText((String) str.subSequence(0, idx));
                                }
                                printOUT(str);
                            }
                        } else {
                            //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
                            printERR("Error: " + retmsg.getRetMSG() + " !");
                        }
                    }
                } else { // Le choix n'a pas été fait
                    printERR("Choose a destination option please at config. tab !");
                }
                abinitPathButton.setEnabled(true);
            }
        };

        Thread t = new Thread(r);
        t.start();
}//GEN-LAST:event_abinitPathButtonActionPerformed

    private void parallelCheckBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_parallelCheckBoxActionPerformed
        parallelCheckBox.setForeground(Color.red);
        sequentialCheckBox.setForeground(Color.blue);

        abinitParaTextField.setEnabled(true);
        abinitParaLabel.setEnabled(true);
}//GEN-LAST:event_parallelCheckBoxActionPerformed

    private void sequentialCheckBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_sequentialCheckBoxActionPerformed
        sequentialCheckBox.setForeground(Color.red);
        parallelCheckBox.setForeground(Color.blue);

        abinitParaTextField.setEnabled(false);
        abinitParaLabel.setEnabled(false);
}//GEN-LAST:event_sequentialCheckBoxActionPerformed

    private void connectionToggleButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_connectionToggleButtonActionPerformed
        Runnable r = new Runnable() {

            @Override
            public void run() {
                connectionToggleButton.setEnabled(false);
                if (connectionToggleButton.isSelected()) {
                    if (remoteGatewayRadioButton.isSelected()) {
                        String gwHostname = gatewayHostTextField.getText();
                        String gwLogin = gatewayLoginTextField.getText();
                        if (!gwLogin.equals("") && !gwHostname.equals("")) {
                            String abHostname = hostTextField.getText();
                            String abLogin = loginTextField.getText();
                            if (!abLogin.equals("") && !abHostname.equals("")) {
                                // Début de la création du tunnel SSH
                                printOUT("Connecting to " + gwHostname + " as " + gwLogin + ".");
                                sshtun = new SSHTunnel(gwLogin, gwHostname, 22, abHostname, 2568, 22);
                                sshtun.setDialog(outDialog);
                                sshtun.setPassword(new String(gatewayPasswordField.getPassword()));
                                lport = sshtun.start();
                                if (lport > 0 && lport < 65536) {
                                    printOUT("Connected to " + gwHostname + " as " + gwLogin + ".");
                                    printOUT("Connecting to " + abHostname + " as " + abLogin + ".");
                                    remoteExec = new RemoteExec(abLogin, "localhost", lport);
                                    remoteExec.setDialog(outDialog);
                                    remoteExec.setPassword(new String(pwdPasswordField.getPassword()));
                                    if (remoteExec.start()) {
                                        printOUT("Connected to " + abHostname + " as " + abLogin + ".");
                                        // Le tunnel SSH a été créé avec succÃ¨s
                                        connectionToggleButton.setText("Disconnect");
                                    } else {
                                        lport = 0;
                                        printERR("Could not connect to " + abHostname + " as " + abLogin + " !");
                                        connectionToggleButton.setSelected(false);
                                        stopConnection();
                                    }
                                } else {
                                    lport = 0;
                                    printERR("Could not connect to " + gwHostname + " as " + gwLogin + " !");
                                    connectionToggleButton.setSelected(false);
                                    stopConnection();
                                }
                            } else {
                                printERR("Please enter the ABINIT hostname AND corresponding login !");
                                connectionToggleButton.setSelected(false);
                                stopConnection();
                            }
                        } else {
                            printERR("Please enter the gateway hostname AND corresponding login !");
                            connectionToggleButton.setSelected(false);
                            stopConnection();
                        }
                    } else if (remoteAbinitRadioButton.isSelected()) {
                        String abHostname = hostTextField.getText();
                        String abLogin = loginTextField.getText();
                        if (!abLogin.equals("") && !abHostname.equals("")) {
                            printOUT("Connecting to " + abHostname + " as " + abLogin + ".");
                            remoteExec = new RemoteExec(abLogin, abHostname, 22);
                            remoteExec.setDialog(outDialog);
                            remoteExec.setPassword(new String(pwdPasswordField.getPassword()));
                            if (remoteExec.start()) {
                                printOUT("Connected to " + abHostname + " as " + abLogin + ".");
                                connectionToggleButton.setText("Disconnect");
                            } else {
                                printERR("Could not connect to " + abHostname + " as " + abLogin + " !");
                                connectionToggleButton.setSelected(false);
                                stopConnection();
                            }
                        } else {
                            printERR("Please enter the ABINIT hostname AND corresponding login !");
                            connectionToggleButton.setSelected(false);
                            stopConnection();
                        }
                    } else if (localAbinitRadioButton.isSelected()) {
                        // Pas besoin en local
                    } else { // Le choix n'a pas été fait
                        printERR("Choose a destination option please at config. tab !");
                    }
                } else {
                    stopConnection();
                    connectionToggleButton.setText("Connect");
                    printOUT("You are now unconnected!");
                }
                connectionToggleButton.setEnabled(true);
            }
        };

        Thread t = new Thread(r);
        t.start();
}//GEN-LAST:event_connectionToggleButtonActionPerformed

    private void needSGECheckBoxActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_needSGECheckBoxActionPerformed
        // TODO add your handling code here:
        if (needSGECheckBox.isSelected()) {
            SGEconfigPanel.setVisible(true);
        } else {
            SGEconfigPanel.setVisible(false);
        }
}//GEN-LAST:event_needSGECheckBoxActionPerformed

    private void remoteGatewayRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_remoteGatewayRadioButtonActionPerformed
        connectionToggleButton.setEnabled(true);

        localAbinitRadioButton.setForeground(Color.blue);
        remoteGatewayRadioButton.setForeground(Color.red);
        remoteAbinitRadioButton.setForeground(Color.blue);

        gatewayLoginPanel.setVisible(true);
        loginPanel.setVisible(true);
        loginPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(
                javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)),
                "Remote Abinithost login", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION,
                javax.swing.border.TitledBorder.DEFAULT_POSITION,
                new java.awt.Font("Arial", 3, 14), java.awt.Color.darkGray));
}//GEN-LAST:event_remoteGatewayRadioButtonActionPerformed

    private void remoteAbinitRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_remoteAbinitRadioButtonActionPerformed
        connectionToggleButton.setEnabled(true);

        localAbinitRadioButton.setForeground(Color.blue);
        remoteGatewayRadioButton.setForeground(Color.blue);
        remoteAbinitRadioButton.setForeground(Color.red);

        gatewayLoginPanel.setVisible(false);
        loginPanel.setVisible(true);
        loginPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(
                javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)),
                "Remote Abinithost login", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION,
                javax.swing.border.TitledBorder.DEFAULT_POSITION,
                new java.awt.Font("Arial", 3, 14), java.awt.Color.darkGray));
}//GEN-LAST:event_remoteAbinitRadioButtonActionPerformed

    private void localAbinitRadioButtonActionPerformed(java.awt.event.ActionEvent evt) {//GEN-FIRST:event_localAbinitRadioButtonActionPerformed
        connectionToggleButton.setEnabled(false);

        localAbinitRadioButton.setForeground(Color.red);
        remoteGatewayRadioButton.setForeground(Color.blue);
        remoteAbinitRadioButton.setForeground(Color.blue);

        gatewayLoginPanel.setVisible(false);
        loginPanel.setVisible(false);
        loginPanel.setBorder(javax.swing.BorderFactory.createTitledBorder(
                javax.swing.BorderFactory.createLineBorder(new java.awt.Color(0, 0, 0)),
                "Local Abinithost login", javax.swing.border.TitledBorder.DEFAULT_JUSTIFICATION,
                javax.swing.border.TitledBorder.DEFAULT_POSITION,
                new java.awt.Font("Arial", 3, 14), java.awt.Color.darkGray));
}//GEN-LAST:event_localAbinitRadioButtonActionPerformed

    private void mySimulationsLabelMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_mySimulationsLabelMouseClicked
        printGEN("--- HINT ------------------------------------------", Color.BLACK, false, true);
        printGEN("You have to start your path with ./ and give a folder name where"
                + " to create the abinit filetree\n", Color.RED, false, true);
        printGEN("Example: ./MySimulations\n", new Color(0,100,0), false, true);
        printGEN("The filetree will be created in your local computer and at the"
                + " Abinit server side when using remote Abinit servers", Color.DARK_GRAY, false, true);
        printGEN("---------------------------------------------------", Color.BLACK, false, true);
    }//GEN-LAST:event_mySimulationsLabelMouseClicked

    private void pspPathLabelMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_pspPathLabelMouseClicked
        printGEN("--- HINT ------------------------------------------", Color.BLACK, false, true);
        printGEN("Please fill in the path to your local pseudopotential database\n", Color.RED, false, true);
        printGEN("Examples: ./PSP (when the database is in your root folder or in the"
                + " same one as AbinitGUI.jar) or something like /home/user/PSP\n", new Color(0,100,0), false, true);
        printGEN("You can find the database at http://www.flavio-abreu.net", Color.DARK_GRAY, false, true);
        printGEN("---------------------------------------------------", Color.BLACK, false, true);
    }//GEN-LAST:event_pspPathLabelMouseClicked

    private void abinitPathPathLabelMouseClicked(java.awt.event.MouseEvent evt) {//GEN-FIRST:event_abinitPathPathLabelMouseClicked
        printGEN("--- HINT ------------------------------------------", Color.BLACK, false, true);
        printGEN("Remote path where to find the abinit program\n", Color.RED, false, true);
        printGEN("Example: /Users/me/Abinit6.7.2/bin\n", new Color(0,100,0), false, true);
        printGEN("---------------------------------------------------", Color.BLACK, false, true);
    }//GEN-LAST:event_abinitPathPathLabelMouseClicked

    private void sendCommand(String CMD) /*throws CMDException*/ {
        RetMSG retmsg;
        if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
            retmsg = remoteExec.sendCommand(CMD);
            if (retmsg.getRetCode() == RetMSG.SUCCES) {
                printOUT("Succes: " + retmsg.getCMD() + " => " + removeEndl(retmsg.getRetMSG()) + ".");
            } else {
                //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
                printERR("Error: " + removeEndl(retmsg.getRetMSG()) + " !");
            }
        } else if (localAbinitRadioButton.isSelected()) {
            retmsg = localExec.sendCommand(CMD);
            if (retmsg.getRetCode() == RetMSG.SUCCES) {
                printOUT("Succes: " + retmsg.getCMD() + " => " + removeEndl(retmsg.getRetMSG()) + ".");
            } else {
                //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
                printERR("Error: " + removeEndl(retmsg.getRetMSG()) + " !");
            }
        } else { // Le choix n'a pas été fait
            printERR("Choose a destination option please at config. tab !");
        }
    }

    private void localCommand(String CMD) /*throws CMDException*/ {
        RetMSG retmsg;
        retmsg = localExec.sendCommand(CMD);
        if (retmsg.getRetCode() == RetMSG.SUCCES) {
            printOUT("Succes: " + retmsg.getCMD() + " => " + removeEndl(retmsg.getRetMSG()) + ".");
        } else {
            //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
            printERR("Error: " + removeEndl(retmsg.getRetMSG()) + " !");
        }
    }

    private void mkdir(String dir) {
        if (Utils.mkdir(dir)) {
            printOUT("Succes: mkdir " + dir + ".");
        } else {
            if (Utils.exists(dir)) {
                printDEB("The local directory `" + dir + "' exists !");
            } else {
                printERR("Error: mkdir: cannot create directory `" + dir + "' !");
            }
        }
    }

    private void mkdirR(String dir) {
        String CMD = "mkdir " + dir;
        RetMSG retmsg;
        if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
            retmsg = remoteExec.sendCommand(CMD);
            if (retmsg.getRetCode() == RetMSG.SUCCES) {
                printOUT("Succes: " + retmsg.getCMD() + " => " + removeEndl(retmsg.getRetMSG()) + ".");
            } else {
                if (retmsg.getRetCode() == 1) {
                    printDEB("The remote directory `" + dir + "' exists !");
                } else {
                    //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
                    printERR("Error: " + removeEndl(retmsg.getRetMSG()) + " !");
                }
            }
        } else {
            printERR("Error: local use of mkdirR !");
        }
    }

    private String getOutputFilesR(String dir) {
        String CMD = "ls " + dir;
        RetMSG retmsg;
        if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
            retmsg = remoteExec.sendCommand(CMD);
            if (retmsg.getRetCode() == RetMSG.SUCCES) {
                printOUT("Succes: " + retmsg.getCMD());
                return removeEndl(retmsg.getRetMSG());
            } else {
                //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
                printERR("Error: " + removeEndl(retmsg.getRetMSG()) + " !");
                return "";
            }
        } else {
            printERR("Error: local use of lsR !");
            return "";
        }
    }

    private void localCopy(String parameters) {
        RetMSG retmsg = localExec.sendCommand("cp " + parameters);
        System.out.println(retmsg.getRetMSG());
        if (retmsg.getRetCode() == RetMSG.SUCCES) {
            printOUT("Succes: " + retmsg.getCMD() + " => " + removeEndl(retmsg.getRetMSG()) + ".");
        } else {
            //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
            printERR("Error: " + removeEndl(retmsg.getRetMSG()) + " !");
        }
    }

    private void putFile(String parameters) {
        RetMSG retmsg;
        if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
            retmsg = remoteExec.sendCommand("put " + parameters);
            if (retmsg.getRetCode() == RetMSG.SUCCES) {
                printOUT("Succes: " + retmsg.getCMD() + " => " + removeEndl(retmsg.getRetMSG()) + ".");
            } else {
                //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
                printERR("Error: " + removeEndl(retmsg.getRetMSG()) + " !");
            }
        } else if (localAbinitRadioButton.isSelected()) {
            retmsg = localExec.sendCommand("cp " + parameters);
            if (retmsg.getRetCode() == RetMSG.SUCCES) {
                printOUT("Succes: " + retmsg.getCMD() + " => " + removeEndl(retmsg.getRetMSG()) + ".");
            } else {
                //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
                printERR("Error: " + removeEndl(retmsg.getRetMSG()) + " !");
            }
        } else { // Le choix n'a pas été fait
            printERR("Choose a destination option please at config. tab !");
        }
    }

    private void getFile(String parameters) {
        RetMSG retmsg;
        if (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected()) {
            retmsg = remoteExec.sendCommand("get " + parameters);
            if (retmsg.getRetCode() == RetMSG.SUCCES) {
                printOUT("Succes: " + retmsg.getCMD() + " => " + removeEndl(retmsg.getRetMSG()) + ".");
            } else {
                //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
                printERR("Error: " + removeEndl(retmsg.getRetMSG()) + " !");
            }
        } else if (localAbinitRadioButton.isSelected()) {
            retmsg = localExec.sendCommand("cp " + parameters);
            if (retmsg.getRetCode() == RetMSG.SUCCES) {
                printOUT("Succes: " + retmsg.getCMD() + " => " + removeEndl(retmsg.getRetMSG()) + ".");
            } else {
                //printERR("Error (RetVal = " + retmsg.getRetCode() + "): " + retmsg.getRetMSG());
                printERR("Error: " + removeEndl(retmsg.getRetMSG()) + " !");
            }
        } else { // Le choix n'a pas été fait
            printERR("Choose a destination option please at config. tab !");
        }
    }

    public class SFTPCMDThread extends Thread {

        private String CMD;

        public SFTPCMDThread(String CMD) {
            this.CMD = CMD;
            start();
        }

        @Override
        public void run() {
            sftp.sendCommand(CMD);
            SFTPCommandLine.setText("");
            sendSFTPButton.setEnabled(true);
            SFTPCommandLine.setEditable(true);
        }
    }

    public class SendSFTPThread extends Thread {

        private String path_;

        public SendSFTPThread(String path) {
            this.path_ = path;
            start();
        }

        @Override
        public void run() {
            JFileChooser fc = new JFileChooser(path_);
            fc.setMultiSelectionEnabled(false);

            int retValue = fc.showDialog(null, "Send");//.showOpenDialog(this);

            if (retValue == JFileChooser.APPROVE_OPTION) {
                File file = fc.getSelectedFile();
                String path = file.getPath();
                try {
                    sftp.sendCommand("put " + path);
                } catch (Exception e) {
                    printERR(e.getMessage());
                }
            }
            sendSFTPButton.setEnabled(true);
            SFTPCommandLine.setEditable(true);
        }
    }

    void saveConfig(String confFile) {
        String tmpvar = "";

        XMLConfigWriter conf = new XMLConfigWriter("tabbedpaneconf");

        // ********************* congfiguration ********************************
        // *********************************************************************
        Element congfiguration = conf.add2root("configuration");
        //**********************************************************************
        Enumeration en = whereIsAbinitbuttonGroup.getElements();
        while (en.hasMoreElements()) {
            JRadioButton jrb = (JRadioButton) en.nextElement();
            if (jrb.isSelected()) {
                tmpvar = jrb.getText();
            }
        }
        conf.setAttr(congfiguration, "hostlocation", "lacation", tmpvar);
        //**********************************************************************
        if (needSGECheckBox.isSelected()) {
            tmpvar = "1";
        } else {
            tmpvar = "0";
        }
        conf.setAttr(congfiguration, "needpbs", "checked", tmpvar);
        //**********************************************************************
        tmpvar = loginTextField.getText();
        conf.setAttr(congfiguration, "remotelogin", "login", tmpvar);
        //**********************************************************************
        tmpvar = hostTextField.getText();
        conf.setAttr(congfiguration, "remotehost", "host", tmpvar);
        //**********************************************************************
        tmpvar = new String(pwdPasswordField.getPassword());
        conf.setAttr(congfiguration, "remotepwd", "pwd", tmpvar);
        //**********************************************************************
        tmpvar = gatewayLoginTextField.getText();
        conf.setAttr(congfiguration, "gatewaylogin", "login", tmpvar);
        //**********************************************************************
        tmpvar = gatewayHostTextField.getText();
        conf.setAttr(congfiguration, "gatewayhost", "host", tmpvar);
        //**********************************************************************
        tmpvar = new String(gatewayPasswordField.getPassword());
        conf.setAttr(congfiguration, "gatewaypwd", "pwd", tmpvar);
        //**********************************************************************
        tmpvar = mySimulationsTextField.getText();
        conf.setAttr(congfiguration, "mySimulationspath", "path", tmpvar);
        //**********************************************************************
        tmpvar = pspPathTextField.getText();
        conf.setAttr(congfiguration, "psppath", "path", tmpvar);
        //**********************************************************************
        tmpvar = abinitPathTextField.getText();
        conf.setAttr(congfiguration, "abinitpath", "path", tmpvar);
        //**********************************************************************
        if (sequentialCheckBox.isSelected()) {
            tmpvar = "1";
        } else {
            tmpvar = "0";
        }
        conf.setAttr(congfiguration, "abinitSequ", "checked", tmpvar);
        //**********************************************************************
        if (parallelCheckBox.isSelected()) {
            tmpvar = "1";
        } else {
            tmpvar = "0";
        }
        conf.setAttr(congfiguration, "abinitPara", "checked", tmpvar);
        //**********************************************************************
        tmpvar = abinitParaTextField.getText();
        conf.setAttr(congfiguration, "abinipproc", "numproc", tmpvar);
        //**********************************************************************
        tmpvar = timeTextField.getText();
        conf.setAttr(congfiguration, "pbstime", "time", tmpvar);
        //**********************************************************************
        tmpvar = nodesTextField.getText();
        conf.setAttr(congfiguration, "pbsnodes", "nodes", tmpvar);
        //**********************************************************************
        tmpvar = ramTextField.getText();
        conf.setAttr(congfiguration, "pbsram", "ram", tmpvar);
        //**********************************************************************
        tmpvar = hdmTextField.getText();
        conf.setAttr(congfiguration, "pbshdm", "hdm", tmpvar);
        //**********************************************************************
        tmpvar = emailTextField.getText();
        conf.setAttr(congfiguration, "pbsemail", "email", tmpvar);
        //**********************************************************************
        tmpvar = openFileTextField.getText();
        conf.setAttr(congfiguration, "openfile", "file", tmpvar);
        //**********************************************************************
        tmpvar = pspTextField.getText();
        conf.setAttr(congfiguration, "psp", "psp", tmpvar);
        //**********************************************************************
        int nbelem = 0;
        try {
            nbelem = Integer.parseInt(pspTextField.getText());
        } catch (Exception e) {
            nbelem = 0;
        }
        if (nbelem > 0) {
            Element pspTable_ = new Element("psptable");
            congfiguration.addContent(pspTable_);

            int row = pspTable.getRowCount();
            if (row > 0) {
                for (int i = 0; i < row; i++) {

                    try {
                        // L'élément pspTable.getValueAt(i, 0) est de type Atom
                        // et toutes les infos du pseudopotentiel sont dedans.
                        // J'aurais pu utiliser cet objet pour obtenir tous les champs!
                        String symbol = (pspTable.getValueAt(i, 0)).toString();
                        String pspfile = (pspTable.getValueAt(i, 1)).toString();
                        String psptype = (pspTable.getValueAt(i, 2)).toString();
                        String psppath = (pspTable.getValueAt(i, 3)).toString();

                        conf.setAttr(pspTable_, "entry", new String[]{"i", "symbol", "psppath", "pspfile", "psptype"},
                                new String[]{Integer.toString(i), symbol, psppath, pspfile, psptype});
                    } catch (Exception e) {
                        printERR("Bug reading pspTable (XML) !");
                    }

                }
            }
        }

        // ********************* input file ************************************
        // *********************************************************************
        Element inputfile = conf.add2root("inputfile");

        // ********************* ssh terminal **********************************
        // *********************************************************************
        Element sshterminal = conf.add2root("sshterminal");
        //**********************************************************************
        tmpvar = SSHUserAndHostTextField.getText();
        conf.setAttr(sshterminal, "sshuserhost", "userhost", tmpvar);
        //**********************************************************************
        tmpvar = new String(SSHPwdPasswordField.getPassword());
        conf.setAttr(sshterminal, "sshpwd", "pwd", tmpvar);

        // ********************* sftp terminal *********************************
        // *********************************************************************
        Element sftpterminal = conf.add2root("sftpterminal");
        //**********************************************************************
        tmpvar = SFTPUserAndHostTextField.getText();
        conf.setAttr(sftpterminal, "sftpuserhost", "userhost", tmpvar);
        //**********************************************************************
        tmpvar = new String(SFTPPwdPasswordField.getPassword());
        conf.setAttr(sftpterminal, "sftppwd", "pwd", tmpvar);

        // *********************************************************************
        // Pour le débuggage
        //conf.display();

        if (confFile == null) {
            conf.save2file("config.xml");
        } else {
            conf.save2file(confFile);
        }
    }

    void loadConfig(String file2load) {

        XMLConfigReader conf = new XMLConfigReader(file2load);

        if (conf.getRoot() != null) {
            List l1 = conf.getRoot().getChildren();
            Iterator i1 = l1.iterator();
            while (i1.hasNext()) {
                Element cur1 = (Element) i1.next();
                // *************** configuration ***********************************
                if (cur1.getName().equals("configuration")) {
                    List l2 = cur1.getChildren();
                    Iterator i2 = l2.iterator();
                    while (i2.hasNext()) {
                        Element cur2 = (Element) i2.next();
                        String elemName = cur2.getName();
                        Attribute attr;
                        List lAttr = cur2.getAttributes();
                        Iterator iAttr = lAttr.iterator();
                        // On s'intéresse qu'Ã  un seul atribut
                        if (iAttr.hasNext()) {
                            attr = (Attribute) iAttr.next();
                            String attrValue = attr.getValue();

                            if (elemName.equals("hostlocation")) {
                                if (attrValue.equals(localAbinitRadioButton.getText())) {
                                    localAbinitRadioButton.setSelected(true);
                                    localAbinitRadioButtonActionPerformed(null);
                                } else if (attrValue.equals(remoteAbinitRadioButton.getText())) {
                                    remoteAbinitRadioButton.setSelected(true);
                                    remoteAbinitRadioButtonActionPerformed(null);
                                } else if (attrValue.equals(remoteGatewayRadioButton.getText())) {
                                    remoteGatewayRadioButton.setSelected(true);
                                    remoteGatewayRadioButtonActionPerformed(null);
                                }
                            } else if (elemName.equals("needpbs")) {
                                if (attrValue.equals("1")) {
                                    needSGECheckBox.setSelected(true);
                                    SGEconfigPanel.setVisible(true);
                                } else {
                                    needSGECheckBox.setSelected(false);
                                    SGEconfigPanel.setVisible(false);
                                }
                            } else if (elemName.equals("remotelogin")) {
                                loginTextField.setText(attrValue);
                            } else if (elemName.equals("remotehost")) {
                                hostTextField.setText(attrValue);
                            } else if (elemName.equals("remotepwd")) {
                                pwdPasswordField.setText(attrValue);
                            } else if (elemName.equals("gatewaylogin")) {
                                gatewayLoginTextField.setText(attrValue);
                            } else if (elemName.equals("gatewayhost")) {
                                gatewayHostTextField.setText(attrValue);
                            } else if (elemName.equals("gatewaypwd")) {
                                gatewayPasswordField.setText(attrValue);
                            } else if (elemName.equals("mySimulationspath")) {
                                mySimulationsTextField.setText(attrValue);
                            } else if (elemName.equals("psppath")) {
                                pspPathTextField.setText(attrValue);
                            } else if (elemName.equals("abinitpath")) {
                                abinitPathTextField.setText(attrValue);
                            } else if (elemName.equals("abinitSequ")) {
                                if (attrValue.equals("1")) {
                                    sequentialCheckBox.setSelected(true);
                                    sequentialCheckBoxActionPerformed(null);
                                }
                            } else if (elemName.equals("abinitPara")) {
                                if (attrValue.equals("1")) {
                                    parallelCheckBox.setSelected(true);
                                    parallelCheckBoxActionPerformed(null);
                                }
                            } else if (elemName.equals("abinipproc")) {
                                abinitParaTextField.setText(attrValue);
                            } else if (elemName.equals("pbstime")) {
                                timeTextField.setText(attrValue);
                            } else if (elemName.equals("pbsnodes")) {
                                nodesTextField.setText(attrValue);
                            } else if (elemName.equals("pbsram")) {
                                ramTextField.setText(attrValue);
                            } else if (elemName.equals("pbshdm")) {
                                hdmTextField.setText(attrValue);
                            } else if (elemName.equals("pbsemail")) {
                                emailTextField.setText(attrValue);
                            } else if (elemName.equals("openfile")) {
                                openFileTextField.setText(attrValue);
                            } else if (elemName.equals("psp")) {
                                pspTextField.setText(attrValue);
                            } else {
                                printERR("Unknown element called " + elemName);
                            }
                        } else {
                            if (elemName.equals("psptable")) {
                                try {
                                    int npsp = Integer.parseInt(pspTextField.getText());

                                    if (npsp > 1000) {
                                        npsp = 1000;
                                        Object strTab[][] = new Object[npsp][4];
                                        for (int i = 0; i < npsp; i++) {
                                            strTab[i] = new Object[]{new Atom(), "", "", ""};
                                        }
                                        pspModel.setData(strTab);
                                        //znuclTable.setModel(znuclModel);
                                    } else {
                                        Object strTab[][] = new Object[npsp][4];
                                        for (int i = 0; i < npsp; i++) {
                                            strTab[i] = new Object[]{new Atom(), "", "", ""};
                                        }
                                        pspModel.setData(strTab);
                                        //znuclTable.setModel(znuclModel);
                                    }
                                } catch (Exception e) {
                                    //printERR(e.getMessage());
                                    pspModel.setData(null);
                                    //znuclTable.setModel(znuclModel);
                                }
                                List l3 = cur2.getChildren();
                                Iterator i3 = l3.iterator();
                                while (i3.hasNext()) {
                                    Element cur3 = (Element) i3.next();
                                    String elemName2 = cur3.getName();

                                    if (elemName2.equals("entry")) {
                                        int i = -1;
                                        String symbol = "";
                                        String psppath = "";
                                        String pspfile = "";
                                        String psptype = "";

                                        Attribute attr2;
                                        List lAttr2 = cur3.getAttributes();
                                        Iterator iAttr2 = lAttr2.iterator();

                                        while (iAttr2.hasNext()) {
                                            attr2 = (Attribute) iAttr2.next();
                                            String attrName = attr2.getName();
                                            String attrValue = attr2.getValue();

                                            if (attrName.equals("i")) {
                                                // TODO géréer une évantuelle exception
                                                i = Integer.parseInt(attrValue);
                                            } else if (attrName.equals("symbol")) {
                                                symbol = attrValue;
                                            } else if (attrName.equals("psppath")) {
                                                psppath = attrValue;
                                            } else if (attrName.equals("pspfile")) {
                                                pspfile = attrValue;
                                            } else if (attrName.equals("psptype")) {
                                                psptype = attrValue;
                                            } else {
                                                printERR("Unknown attribute called " + attrName);
                                            }
                                        }

                                        if (i != -1) {
                                            printDEB(i + " " + symbol + " " + psppath + " " + pspfile + " " + psptype);

                                            Atom atom = (Atom) pspTable.getValueAt(i, 0);

                                            AtomEditor.setAtom(atom, symbol, psptype, psppath, pspfile);
                                        } else {
                                            printERR("The xml entry (psptable) is not well defined !");
                                        }
                                    } else {
                                        printERR("Unknown element called " + elemName2);
                                    }
                                }
                            } else {
                                printERR("There is no attribute for " + elemName);
                            }
                        }
                    }
                } else // *************** inputfile ***************************************
                if (cur1.getName().equals("inputfile")) {
                    //printOUT(cur1.getName());
                } else // *************** ssh terminal ************************************
                if (cur1.getName().equals("sshterminal")) {
                    List l2 = cur1.getChildren();
                    Iterator i2 = l2.iterator();
                    while (i2.hasNext()) {
                        Element cur2 = (Element) i2.next();
                        Attribute attr;
                        List lAttr = cur2.getAttributes();
                        Iterator iAttr = lAttr.iterator();
                        // On s'intéresse qu'Ã  un seul atribut
                        if (iAttr.hasNext()) {
                            attr = (Attribute) iAttr.next();

                            String elemName = cur2.getName();
                            String attrValue = attr.getValue();

                            if (elemName.equals("sshuserhost")) {
                                SSHUserAndHostTextField.setText(attrValue);
                            } else if (elemName.equals("sshpwd")) {
                                SSHPwdPasswordField.setText(attrValue);
                            } else {
                                printERR("Unknown element called " + elemName);
                            }
                        } else {
                            printERR("There is no attribute for " + cur2.getName());
                        }
                    }
                } else // *************** sftp terminal *******************************
                if (cur1.getName().equals("sftpterminal")) {
                    List l2 = cur1.getChildren();
                    Iterator i2 = l2.iterator();
                    while (i2.hasNext()) {
                        Element cur2 = (Element) i2.next();
                        Attribute attr;
                        List lAttr = cur2.getAttributes();
                        Iterator iAttr = lAttr.iterator();
                        // On s'intéresse qu'Ã  un seul atribut
                        if (iAttr.hasNext()) {
                            attr = (Attribute) iAttr.next();

                            String elemName = cur2.getName();
                            String attrValue = attr.getValue();

                            if (elemName.equals("sftpuserhost")) {
                                SFTPUserAndHostTextField.setText(attrValue);
                            } else if (elemName.equals("sftppwd")) {
                                SFTPPwdPasswordField.setText(attrValue);
                            } else {
                                printERR("Unknown element called " + elemName);
                            }
                        } else {
                            printERR("There is no attribute for " + cur2.getName());
                        }
                    }
                } else {
                    printERR("Unknown configuration section " + cur1.getName());
                }
            }
        }
    }

    private void createFiletree() {
        if (localExec != null) {
            String path = mySimulationsTextField.getText();
            if (path.equals("")) {
                path = ".";
            }
            if (localAbinitRadioButton.isSelected()) {
                // Création de l'arborescence locale
                mkdir(path);
                mkdir(path + "/input");
                mkdir(path + "/output");
                mkdir(path + "/wholedata");
                mkdir(path + "/logfiles");
                mkdir(path + "/pseudopot");
            } else {
                if (remoteExec != null && (remoteGatewayRadioButton.isSelected() || remoteAbinitRadioButton.isSelected())) {
                    mkdirR(path);
                    mkdirR(path + "/input");
                    mkdirR(path + "/output");
                    mkdirR(path + "/wholedata");
                    mkdirR(path + "/logfiles");
                    mkdirR(path + "/pseudopot");
                    // Création de l'arborescence locale
                    mkdir(path);
                    mkdir(path + "/input");
                    mkdir(path + "/output");
                    mkdir(path + "/wholedata");
                    mkdir(path + "/logfiles");
                    mkdir(path + "/pseudopot");
                } else {
                    printERR("Please connect to an abinit Host before creating the remote filetree !");
                }
            }
        } else {
            printERR("Program bug in createFiletreeButtonActionPerformed (localExec == null) !");
        }
    }
    // Variables declaration - do not modify//GEN-BEGIN:variables
    javax.swing.JMenuItem LoadMenuItem;
    javax.swing.JTextField SFTPCommandLine;
    javax.swing.JLabel SFTPCommandLineLabel;
    javax.swing.JPanel SFTPLoginPanel;
    javax.swing.JPanel SFTPLoginPanel1;
    javax.swing.JTextArea SFTPOutput;
    javax.swing.JLabel SFTPOutputLabel;
    javax.swing.JScrollPane SFTPOutputScrollPane;
    javax.swing.JProgressBar SFTPProgressBar;
    javax.swing.JLabel SFTPPwdLabel;
    javax.swing.JPasswordField SFTPPwdPasswordField;
    javax.swing.JLabel SFTPUserAndHostLabel;
    javax.swing.JTextField SFTPUserAndHostTextField;
    javax.swing.JPanel SGEconfigPanel;
    javax.swing.JTextField SSHCommandLine;
    javax.swing.JLabel SSHCommandLineLabel;
    javax.swing.JTextArea SSHOutput;
    javax.swing.JLabel SSHOutputLabel;
    javax.swing.JScrollPane SSHOutputScrollPane;
    javax.swing.JLabel SSHPwdLabel;
    javax.swing.JPasswordField SSHPwdPasswordField;
    javax.swing.JLabel SSHUserAndHostLabel;
    javax.swing.JTextField SSHUserAndHostTextField;
    javax.swing.JLabel UploadDownloadINFOLabel;
    javax.swing.JLabel UploadDownloadRATELabel;
    javax.swing.JLabel abinitParaLabel;
    javax.swing.JTextField abinitParaTextField;
    javax.swing.JButton abinitPathButton;
    javax.swing.JLabel abinitPathPathLabel;
    javax.swing.JTextField abinitPathTextField;
    javax.swing.ButtonGroup abinixbuttonGroup;
    javax.swing.JMenuItem aboutMenuItem;
    javax.swing.JButton algoAndConvButton;
    javax.swing.JPanel basicsPanel;
    javax.swing.JScrollPane basicsScrollPane;
    javax.swing.JMenuItem clearOutMSGMenuItem;
    javax.swing.JPanel configPanel;
    javax.swing.JToggleButton connectionToggleButton;
    javax.swing.JButton createButton;
    javax.swing.JButton displayFileButton;
    javax.swing.JLabel emailLabel;
    javax.swing.JTextField emailTextField;
    javax.swing.JPanel emptyPanel;
    javax.swing.JMenu fileMenu;
    javax.swing.JTextField gatewayHostTextField;
    javax.swing.JPanel gatewayLoginPanel;
    javax.swing.JTextField gatewayLoginTextField;
    javax.swing.JPasswordField gatewayPasswordField;
    javax.swing.JButton geditButton;
    javax.swing.JButton geometryButton;
    javax.swing.JMenuItem getLogFileMenuItem;
    javax.swing.JMenuItem getOutputFileMenuItem;
    javax.swing.JLabel hdmLabel;
    javax.swing.JTextField hdmTextField;
    javax.swing.JMenu helpMenu;
    javax.swing.JMenuItem helpMenuItem;
    javax.swing.JLabel hostBFELabel;
    javax.swing.JLabel hostLabel;
    javax.swing.JTextField hostTextField;
    javax.swing.JPanel inputFilePanel;
    javax.swing.JTabbedPane inputFileTabbedPane;
    javax.swing.ButtonGroup inputFilebuttonGroup;
    javax.swing.JButton inputOutputButton;
    javax.swing.JLabel jLabel4;
    javax.swing.JScrollPane jScrollPane5;
    javax.swing.JSeparator jSeparator1;
    javax.swing.JRadioButton localAbinitRadioButton;
    javax.swing.JLabel loginBFELabel;
    javax.swing.JLabel loginLabel;
    javax.swing.JPanel loginPanel;
    javax.swing.JTextField loginTextField;
    javax.swing.ButtonGroup lookAndFeelbuttonGroup;
    javax.swing.JMenuBar mainMenuBar;
    javax.swing.JTabbedPane mainTabbedPane;
    javax.swing.JLabel mySimulationsLabel;
    javax.swing.JTextField mySimulationsTextField;
    javax.swing.JCheckBox needSGECheckBox;
    javax.swing.JLabel nodesLabel;
    javax.swing.JTextField nodesTextField;
    javax.swing.JButton openFileDialogButton;
    javax.swing.JLabel openFileLabel;
    javax.swing.JTextField openFileTextField;
    javax.swing.JButton openXMLFileDialogButton;
    javax.swing.JLabel openXMLFileLabel;
    javax.swing.JTextField openXMLFileTextField;
    javax.swing.JTextArea otherTextArea;
    javax.swing.JMenuItem outputMSGMenuItem;
    javax.swing.JCheckBox parallelCheckBox;
    javax.swing.JMenu postProcMenu;
    javax.swing.JLabel pspLabel;
    javax.swing.JLabel pspPathLabel;
    javax.swing.JTextField pspPathTextField;
    javax.swing.JTable pspTable;
    javax.swing.JScrollPane pspTableScrollPane;
    javax.swing.JTextField pspTextField;
    javax.swing.JLabel pwdBFELabel;
    javax.swing.JLabel pwdLabel;
    javax.swing.JPasswordField pwdPasswordField;
    javax.swing.JLabel ramLabel;
    javax.swing.JTextField ramTextField;
    javax.swing.JButton realAndRecipButton;
    javax.swing.JRadioButton remoteAbinitRadioButton;
    javax.swing.JRadioButton remoteGatewayRadioButton;
    javax.swing.JMenuItem saveAsMenuItem;
    javax.swing.JButton saveFileAsButton;
    javax.swing.JButton saveFileButton;
    javax.swing.JMenuItem saveMenuItem;
    javax.swing.JButton sendAFileButton;
    javax.swing.JButton sendSFTPButton;
    javax.swing.JButton sendSIMButton;
    javax.swing.JButton sendSSHButton;
    javax.swing.JCheckBox sequentialCheckBox;
    javax.swing.JPanel sftpPanel;
    javax.swing.JPanel sshPanel;
    javax.swing.JButton startSFTPButton;
    javax.swing.JButton startSSHButton;
    javax.swing.JButton stopSFTPButton;
    javax.swing.JButton stopSSHButton;
    javax.swing.JButton theoryButton;
    javax.swing.JLabel timeLabel;
    javax.swing.JTextField timeTextField;
    javax.swing.JRadioButton useCreIFRadioButton;
    javax.swing.JRadioButton useExtIFRadioButton;
    javax.swing.JCheckBox useGlobalConfigSFTPCheckBox;
    javax.swing.JCheckBox useGlobalConfigSSHCheckBox;
    javax.swing.JMenuItem varsHelpMenuItem;
    javax.swing.JMenu viewMenu;
    javax.swing.JButton wavefuncAndDensButton;
    javax.swing.JLabel whereIsAbinitLabel;
    javax.swing.ButtonGroup whereIsAbinitbuttonGroup;
    // End of variables declaration//GEN-END:variables
}
