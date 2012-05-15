/*--
 SSH.java - Created in July 2009

 Copyright (c) 2009-2011 Flavio Miguel ABREU ARAUJO.
 Universit� catholique de Louvain, Louvain-la-Neuve, Belgium
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

import com.jcraft.jsch.Channel;
import com.jcraft.jsch.JSch;
import com.jcraft.jsch.Session;
import com.jcraft.jsch.UserInfo;
import java.io.InputStream;
import java.io.OutputStream;
import java.util.Vector;
import javax.swing.JTextArea;

//@SuppressWarnings({"deprecation","unchecked"})
public class SSH {

    private InputStream in;
    private OutputStream out;
    public Channel channel;
    public Session session;
    private String userAndHost = null;
    private String password = null;
    private JTextArea outputTA;
    private OutputHandler sshOH;
    private Vector cmds = new Vector();
    private String retMSG = "";
    private SSHExecOutput sshEO;
    boolean console = false;
    private int sshPort = 22;
    private static boolean DEBUG = false;

    private DisplayerJDialog dialog;
    private boolean graphical = false;

    public void setDialog(DisplayerJDialog dialog) {
        graphical = true;
        this.dialog = dialog;
    }

    void printOUT(String str) {
        if (graphical) {
            if (str.endsWith("\n")) {
                dialog.appendOUT(str);
            } else {
                dialog.appendOUT(str + "\n");
            }
        } else {
            if (str.endsWith("\n")) {
                System.out.print(str);
            } else {
                System.out.println(str);
            }
        }
    }

    void printERR(String str) {
        if (graphical) {
            if (str.endsWith("\n")) {
                dialog.appendERR(str);
            } else {
                dialog.appendERR(str + "\n");
            }
        } else {
            if (str.endsWith("\n")) {
                System.err.print(str);
            } else {
                System.err.println(str);
            }
        }
    }

    void printDEB(String str) {
        if (graphical) {
            if (str.endsWith("\n")) {
                dialog.appendDEB(str);
            } else {
                dialog.appendDEB(str + "\n");
            }
        } else {
            if (str.endsWith("\n")) {
                System.out.print(str);
            } else {
                System.out.println(str);
            }
        }
    }

    public SSH(JTextArea outputTA) {
        this(outputTA, false);
    }

    public SSH(JTextArea outputTA, boolean console) {
        this.outputTA = outputTA;
        this.console = console;
    }

    public void setUserAndHost(String userAndHost) {
        this.userAndHost = userAndHost;
    }

    public void setPassword(String password) {
        this.password = password;
    }

    public void setPort(int port) {
        this.sshPort = port;
    }

    /*synchronized*/ public void sendCommand(String cmd) {
        if (console) {
            try {
                cmd += '\n';
                out.write(cmd.getBytes(), 0, cmd.length());
                out.flush();
            } catch (Exception e) {
                printERR(e.getMessage());
            }
        } else {
            if (cmd.startsWith("ls")) {
                cmd += " --color=never";
            }

            if (cmd.startsWith("exec") || cmd.startsWith("EXEC")) {
                // cette commande n'est pas gérée
                outputTA.append("Command [" + cmd + "] not supported by ABINIT GUI !\n");
                return;
            }

            if (cmd.equals("vi") || cmd.equals("VI")) {
                // cette commande n'est pas gérée
                outputTA.append("Command [" + cmd + "] not supported by ABINIT GUI !\n");
                return;
            }

            if (cmd.equals("vim") || cmd.equals("VIM")) {
                // cette commande n'est pas gérée
                outputTA.append("Command [" + cmd + "] not supported by ABINIT GUI !\n");
                return;
            }

            if (cmds.add(cmd)) {
                //printOUT("Command [" + cmd + "] successfuly registred.");
            } else {
                printERR("Could not register the command: " + cmd);
            }
        }
    }

    public boolean start() {
        try {
            if (DEBUG) {
                JSch.setLogger(new MyLogger());
            }
            JSch jsch = new JSch();

            String user = userAndHost.substring(0, userAndHost.indexOf('@'));
            String host = userAndHost.substring(userAndHost.indexOf('@') + 1);

            session = jsch.getSession(user, host, sshPort);

            // username and password will be given via UserInfo interface.
            UserInfo ui = new MyUserInfo();
            ((MyUserInfo) ui).setPassword(password);
            session.setUserInfo(ui);
            session.connect();

            channel = session.openChannel("shell");

            //((ChannelShell) channel).setXForwarding(true); // XForwarding

            out = channel.getOutputStream();
            in = channel.getInputStream();

            if (console) {
                sshEO = new SSHExecOutput(in);
            } else {
                sshOH = new OutputHandler();
                sshOH.start();
            }
            channel.connect();

            if (session.isConnected() && channel.isConnected()) {
                return true;
            } else {
                return false;
            }
        } catch (com.jcraft.jsch.JSchException e) {
            String msg = e.getMessage();
            if (e.getCause() != null) {
                String msg2 = msg.substring(0, msg.indexOf(':'));
                String msg3 = e.getCause().getMessage();
                if (msg2.equals("java.net.UnknownHostException")) {
                    printERR("Unknown hostname: " + msg3);
                } else {
                    printERR("Exception: " + msg2);
                }
            } else {
                if (msg.equals("Auth fail")) {
                    printERR("Username or password is wrong !!");
                } else {
                    printERR("JSchException => MSG: " + msg);
                }
            }
            return false;
        } catch (Exception e) {
            printERR(e.getMessage());
            return false;
        }
    }

    public void stop() {
        if (console) {
            if(channel != null && session != null) sendCommand("exit");
            if (channel != null) {
                channel.disconnect();
            }
            if (session != null) {
                session.disconnect();
            }
            if (sshEO != null) {
                sshEO.stop();
            }
        } else {
            new ClosureWatcher().start();
        }
    }

    public boolean isConnected() {
        if(channel.isConnected()) {
            return session.isConnected();
        } else {
            return false;
        }
    }

    public class OutputHandler extends Thread {

        //private String lastCMD = "";
        @Override
        public void run() {
            String msg = null;
            byte[] buf = new byte[1024];
            int i = 0;
            try {
                while (true) {
                    msg = "";
                    do {
                        int b = in.read();
                        //outputTA.append("[" + b + "-" + (char) b + "]\n");
                        if (b == -1) {
                            //OutputTA.append("code de retour [" + b + "]\n");
                        } else if (b == 0) {
                            outputTA.append("[" + b + "] SUCCES\n");
                        } else if (b == 1) {
                            outputTA.append("[" + b + "] ERROR\n");
                        } else if (b == 2) {
                            outputTA.append("[" + b + "] FATAL ERROR\n");
                        }

                        i = in.read(buf, 0, 1024);
                        if (i <= 0) {
                            break;
                        } else {
                            byte[] tmp = new byte[i];
                            for (int k = 0; k < i; k++) {
                                tmp[k] = buf[k];
                            }

                            msg = (char) b + new String(tmp);
                            outputTA.append(msg);
                            outputTA.setCaretPosition(outputTA.getDocument().getLength());

                            if (msg.contains("@") && ((msg.endsWith("> ") || msg.endsWith(">") ||
                                    msg.endsWith("# ") || msg.endsWith("#") ||
                                    msg.endsWith("$ ") || msg.endsWith("$")))) {
                                // THE PROMPT IS READY !!
                                //printOUT(msg + '\n');
                                String cmd = execNextCMD();
                                if (cmd.startsWith("exit") || cmd.startsWith("EXIT")) {
                                    //Ceci provoque l'arrêt contrôlé du thread OutputHandler
                                    throw new Exception("OutputHandler thread is shuting down !!");
                                }
                            } else {
                                if (msg.contains("Password:") || msg.contains("password:")) {
                                    printOUT("msg.contains(\"Password:\") || msg.contains(\"password:\")");
                                    execNextCMD();
                                } else if (msg.endsWith("ETA")) {
                                    // scp command
                                    outputTA.append("\n");
                                    outputTA.setCaretPosition(outputTA.getDocument().getLength());
                                } else {
                                    // Continuer à lire
                                    break;
                                }
                            }

                        }
                    } while (in.available() > 0);
                }
            } catch (Exception e) {
               printERR(e.getMessage());
            }
        }

        private String execNextCMD() {
            String cmd = null;
            do {
                if (!cmds.isEmpty()) {
                    cmd = (String) cmds.firstElement();
                    break;
                }
                try {
                    Thread.sleep(10);
                } catch (Exception e) {
                }
            } while (true);
            cmds.remove(cmd);
            try {
                cmd += '\n';
                out.write(cmd.getBytes(), 0, cmd.length());
                out.flush();
            } catch (Exception e) {
                printERR(e.getMessage());
            }
            return cmd;
        }
    }

    public class ClosureWatcher extends Thread {

        @Override
        public void run() {
            sendCommand("exit");
            try {
                sshOH.join();
            } catch (Exception e) {
                printERR("SSH stop(): " + e);
            }

            if (channel != null) {
                channel.disconnect();
            }

            if (session != null) {
                session.disconnect();
            }

            if (sshOH != null) {
                sshOH.stop();
            }
        }
    }

    public class SSHExecOutput extends Thread {

        private InputStream input;

        public SSHExecOutput(InputStream input) {
            this.input = input;
            this.start();
        }

        @Override
        public void run() {
            byte[] buf = new byte[1024];
            int i = 0;
            try {
                while (true) {
                    i = input.read(buf, 0, 1024);
                    if (i <= 0) {
                        break;
                    } else {
                        byte[] tmp = new byte[i];
                        for (int k = 0; k < i; k++) {
                            tmp[k] = buf[k];
                        }
                        String str = new String(tmp);

                        outputTA.append(str.replace("[0;0m", "")
                                .replace("[1;39m", "")
                                .replace("[00m", "")
                                .replace("[01;31m", "")
                                .replace("[01;32m", "")
                                .replace("[01;33m", "")
                                .replace("[01;34m", "")
                                .replace("[m", "")
                                .replace("[0m", "")
                                .replace("[K", ""));
                        outputTA.setCaretPosition(outputTA.getDocument().getLength());
                    }
                }
            } catch (Exception e) {
                printERR(e.getMessage());
            }
        }
    }
}


