/*--
 SFTP.java - Created in July 2009

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

import com.jcraft.jsch.Channel;
import com.jcraft.jsch.ChannelSftp;
import com.jcraft.jsch.ChannelSftp.LsEntry;
import com.jcraft.jsch.JSch;
import com.jcraft.jsch.Session;
import com.jcraft.jsch.SftpATTRS;
import com.jcraft.jsch.SftpException;
import com.jcraft.jsch.SftpProgressMonitor;
import com.jcraft.jsch.UserInfo;
import java.util.StringTokenizer;
import javax.swing.JLabel;
import javax.swing.JProgressBar;
import javax.swing.JTextArea;

//@SuppressWarnings("unchecked")
public class SFTP {

    private JProgressBar pb;
    private JLabel lab_file;
    private JLabel lab_rate;
    private JTextArea outputTA;
    private String userAndHost = null;
    private String password = null;
    private ChannelSftp c;
    private Session session;
    private boolean isGraphic;
    private int sftpPort = 22;
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

    void printOUT2(String str) {
        if (graphical) {
            dialog.appendOUT(str);
        } else {
            System.out.print(str);
        }
    }

    public SFTP(JProgressBar pb, JLabel lab_file, JLabel lab_rate, JTextArea outputTA) {
        this.pb = pb;
        this.lab_file = lab_file;
        this.lab_rate = lab_rate;
        this.outputTA = outputTA;
        this.isGraphic = true;
        //this.start();
    }

    /*public SFTP() {
    this(null, null, null, null);
    this.isGraphic = false;
    //this.start();
    }*/
    public void setUserAndHost(String userAndHost) {
        this.userAndHost = userAndHost;
    }

    public void setPassword(String password) {
        this.password = password;
    }

    public boolean isConnected() {
        if (c.isConnected()) {
            return session.isConnected();
        } else {
            return false;
        }
    }

    public void sendCommand(String command) {
        if (!command.equals("")) {
            printText(command + '\n');
        }
        java.util.Vector cmds = new java.util.Vector();
        cmds.removeAllElements();

        StringTokenizer st = new StringTokenizer(command);
        while (st.hasMoreTokens()) {
            cmds.addElement(st.nextToken());
        }

        String str;
        int level = 0;

        if (cmds.size() == 0) {
            return;
        }

        String cmd = (String) cmds.elementAt(0);

        if (cmd.equals("quit")) {
            c.quit();
            // TODO attention !
            return;
        }
        if (cmd.equals("exit")) {
            c.exit();
            // TODO attention !
            return;
        }
        if (cmd.equals("rekey")) {
            try {
                session.rekey();
            } catch (Exception e) {
                printERR(e.getMessage());
            }
            printPrompt();
            return;
        }
        if (cmd.equals("compression")) {
            if (cmds.size() < 2) {
                printText("compression level: " + level + '\n');
                printPrompt();
                return;
            }
            try {
                level = Integer.parseInt((String) cmds.elementAt(1));
                if (level == 0) {
                    session.setConfig("compression.s2c", "none");
                    session.setConfig("compression.c2s", "none");
                } else {
                    session.setConfig("compression.s2c", "zlib@openssh.com,zlib,none");
                    session.setConfig("compression.c2s", "zlib@openssh.com,zlib,none");
                }
                session.rekey();
            } catch (Exception e) {
                printERR(e.getMessage());
            }
            printPrompt();
            return;
        }
        if (cmd.equals("cd") || cmd.equals("lcd")) {
            if (cmds.size() < 2) {
                printPrompt();
                return;
            }
            String path = (String) cmds.elementAt(1);
            if (cmd.equals("cd")) {
                try {
                    c.cd(path);
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            } else {
                try {
                    c.lcd(path);
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            }
            printPrompt();
            return;
        }
        if (cmd.equals("rm") || cmd.equals("rmdir") || cmd.equals("mkdir")) {
            if (cmds.size() < 2) {
                printPrompt();
                return;
            }
            String path = (String) cmds.elementAt(1);
            if (cmd.equals("rm")) {
                try {
                    c.rm(path);
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            } else if (cmd.equals("rmdir")) {
                try {
                    c.rmdir(path);
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            } else {
                try {
                    c.mkdir(path);
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            }
            printPrompt();
            return;
        }
        if (cmd.equals("chgrp") || cmd.equals("chown") || cmd.equals("chmod")) {
            if (cmds.size() != 3) {
                printPrompt();
                return;
            }
            String path = (String) cmds.elementAt(2);
            int foo = 0;
            if (cmd.equals("chmod")) {
                byte[] bar = ((String) cmds.elementAt(1)).getBytes();
                int k;
                for (int j = 0; j < bar.length; j++) {
                    k = bar[j];
                    if (k < '0' || k > '7') {
                        foo = -1;
                        break;
                    }
                    foo <<= 3;
                    foo |= (k - '0');
                }
                if (foo == -1) {
                    printPrompt();
                    return;
                }
            } else {
                try {
                    foo = Integer.parseInt((String) cmds.elementAt(1));
                } catch (Exception e) {
                    printPrompt();
                    return;
                }
            }
            if (cmd.equals("chgrp")) {
                try {
                    c.chgrp(foo, path);
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            } else if (cmd.equals("chown")) {
                try {
                    c.chown(foo, path);
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            } else if (cmd.equals("chmod")) {
                try {
                    c.chmod(foo, path);
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            }
            printPrompt();
            return;
        }

        if (cmd.equals("pwd") || cmd.equals("lpwd")) {
            str = (cmd.equals("pwd") ? "Remote" : "Local");
            str += " working directory: ";
            if (cmd.equals("pwd")) {
                try {
                    str += c.pwd();
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            } else {
                try {
                    str += c.lpwd();
                } catch (Exception e) {
                    printERR(e.getMessage());
                }
            }
            printText(str + '\n');
            printPrompt();
            return;
        }

        if (cmd.equals("ls") || cmd.equals("dir")) {
            String path = ".";
            if (cmds.size() == 2) {
                path = (String) cmds.elementAt(1);
            }
            try {
                java.util.Vector vv = c.ls(path);
                if (vv != null) {
                    for (int ii = 0; ii < vv.size(); ii++) {
                        Object obj = vv.elementAt(ii);
                        if (obj instanceof LsEntry) {
                            printText(((LsEntry) obj).getLongname() + '\n');
                        }
                    }
                }
            } catch (SftpException e) {
                printERR(e.getMessage());
            }
            printPrompt();
            return;
        }

        if (cmd.equals("lls") || cmd.equals("ldir")) {
            String path = ".";
            if (cmds.size() == 2) {
                path = (String) cmds.elementAt(1);
            }
            try {
                java.io.File file = new java.io.File(path);
                if (!file.exists()) {
                    printText(path + ": No such file or directory" + '\n');
                    printPrompt();
                    return;
                }
                if (file.isDirectory()) {
                    String[] list = file.list();
                    for (int ii = 0; ii < list.length; ii++) {
                        printText(list[ii] + '\n');
                    }
                    printPrompt();
                    return;
                }
                printText(path + '\n');
            } catch (Exception e) {
                printERR(e.getMessage());
            }
            printPrompt();
            return;
        }

        if (cmd.equals("get") ||
                cmd.equals("get-resume") || cmd.equals("get-append") ||
                cmd.equals("put") ||
                cmd.equals("put-resume") || cmd.equals("put-append")) {
            if (cmds.size() != 2 && cmds.size() != 3) {
                printERR("Too much parameters for " + (String) cmds.elementAt(0) + " command");
                printPrompt();
                return;
            }
            String p1 = (String) cmds.elementAt(1);
            // String p2=p1;
            String p2 = ".";
            if (cmds.size() == 3) {
                p2 = (String) cmds.elementAt(2);
            }
            SftpProgressMonitor monitor = new MyProgressMonitor(pb, lab_file, lab_rate);
            if (cmd.startsWith("get")) {
                int mode = ChannelSftp.OVERWRITE;
                if (cmd.equals("get-resume")) {
                    mode = ChannelSftp.RESUME;
                } else if (cmd.equals("get-append")) {
                    mode = ChannelSftp.APPEND;
                }
                Get g = new Get(p1, p2, monitor, mode);
                try {
                    g.join();
                } catch (InterruptedException e) {
                    printERR(e.getMessage());
                }
            } else {
                int mode = ChannelSftp.OVERWRITE;
                if (cmd.equals("put-resume")) {
                    mode = ChannelSftp.RESUME;
                } else if (cmd.equals("put-append")) {
                    mode = ChannelSftp.APPEND;
                }
                Put p = new Put(p1, p2, monitor, mode);
                try {
                    p.join();
                } catch (InterruptedException e) {
                    printERR(e.getMessage());
                }
            }
            printPrompt();
            return;
        }

        if (cmd.equals("ln") || cmd.equals("symlink") || cmd.equals("rename")) {
            if (cmds.size() != 3) {
                printPrompt();
                return;
            }
            String p1 = (String) cmds.elementAt(1);
            String p2 = (String) cmds.elementAt(2);
            if (cmd.equals("rename")) {
                try {
                    c.rename(p1, p2);
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            } else {
                try {
                    c.symlink(p1, p2);
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            }
            printPrompt();
            return;
        }

        if (cmd.equals("stat") || cmd.equals("lstat")) {
            if (cmds.size() != 2) {
                printPrompt();
                return;
            }
            String p1 = (String) cmds.elementAt(1);
            SftpATTRS attrs = null;
            if (cmd.equals("stat")) {
                try {
                    attrs = c.stat(p1);
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            } else {
                try {
                    attrs = c.lstat(p1);
                } catch (SftpException e) {
                    printERR(e.getMessage());
                }
            }
            if (attrs != null) {
                printText(attrs.toString() + '\n');
            } else {
            }
            printPrompt();
            return;
        }

        if (cmd.equals("readlink")) {
            if (cmds.size() != 2) {
                printPrompt();
                return;
            }
            String p1 = (String) cmds.elementAt(1);
            String filename = null;
            try {
                filename = c.readlink(p1);
                printText(filename + '\n');
            } catch (SftpException e) {
                printERR(e.getMessage());
            }
            printPrompt();
            return;
        }

        if (cmd.equals("realpath")) {
            if (cmds.size() != 2) {
                printPrompt();
                return;
            }
            String p1 = (String) cmds.elementAt(1);
            String filename = null;
            try {
                filename = c.realpath(p1);
                printText(filename + '\n');
            } catch (SftpException e) {
                printERR(e.getMessage());
            }
            printPrompt();
            return;
        }

        if (cmd.equals("version")) {
            printText("SFTP protocol version " + c.version() + '\n');
            printPrompt();
            return;
        }

        if (cmd.equals("help") || cmd.equals("?")) {
            printText(help + '\n');
            printPrompt();
            return;
        }

        printText("unimplemented command: " + cmd + '\n');
        printPrompt();
    }

    private void printPrompt() {
        if (isGraphic) {
            outputTA.append("sftp> ");
            outputTA.setCaretPosition(outputTA.getDocument().getLength());
        } else {
            printOUT2("sftp> ");
        }
    }

    private void printText(String text) {
        if (isGraphic) {
            outputTA.append(text);
        } else {
            printOUT2(text);
        }
    }

    public ChannelSftp getChannel() {
        return c;
    }

    public void setPort(int port) {
        this.sftpPort = port;
    }

    public boolean start() {
        try {
            if (DEBUG) {
                JSch.setLogger(new MyLogger());
            }
            JSch jsch = new JSch();
            String user = userAndHost.substring(0, userAndHost.indexOf('@'));
            String host = userAndHost.substring(userAndHost.indexOf('@') + 1);
            session = jsch.getSession(user, host, sftpPort);

            // username and password will be given via UserInfo interface.
            UserInfo ui = new MyUserInfo();
            ((MyUserInfo) ui).setPassword(password);
            session.setUserInfo(ui);
            session.connect();
            Channel channel = session.openChannel("sftp");
            channel.connect();
            c = (ChannelSftp) channel;
            printText("sftp> ");
            if (session.isConnected() && c.isConnected()) {
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
        try {
            if (c != null) {
                c.exit();
            }
            if (session != null) {
                session.disconnect();
            }
        } catch (Exception e) {
           printERR("SFTP stop (" + e + ")");
        }
    }

    public class MyProgressMonitor implements SftpProgressMonitor {

        long count = 0;
        long max = 0;
        final long begtime = System.nanoTime();
        long curtime = 0;
        JProgressBar pb_;
        JLabel lab_file_;
        JLabel lab_rate_;

        public MyProgressMonitor(JProgressBar pb, JLabel lab_file, JLabel lab_rate) {
            this.pb_ = pb;
            this.lab_file_ = lab_file;
            this.lab_rate_ = lab_rate;
        }

        @Override
        public void init(int op, String src, String dest, long max) {
            this.max = max;
            lab_file_.setText(((op == SftpProgressMonitor.PUT) ? "PUT" : "GET") +
                    ": " + src + " (" + Utils.inByteFormat(max) + ")");
            count = 0;
            percent = -1;
            pb_.setMaximum((int) max);
            pb_.setValue((int) this.count);
        }
        private long percent = -1;

        @Override
        public boolean count(long count) {
            this.count += count;

            if (percent >= this.count * 100 / max) {
                return true;
            }
            percent = this.count * 100 / max;

            //"Completed " + this.count + "(" + percent + "%) out of " + max + "."
            curtime = System.nanoTime();
            double time = (curtime - begtime) / 1000000000.0;
            double speed = ((double) this.count) / time;

            String timestr;
            if (time < 1.0) {
                timestr = Math.ceil(time * 1000.0) + " ms";
            } else {
                timestr = Utils.formatDouble(time) + " s";
            }

            lab_rate_.setText("Elapsed time " + timestr +
                    " | Speed " + Utils.inByteFormat(speed) + "/s (" + percent + "%)");
            pb_.setValue((int) this.count);

            return true;
        }

        @Override
        public void end() {
            //System.out.println("Transfert completed");
            lab_rate_.setText("<HTML>" + lab_rate_.getText() + " <b>(Transfert completed)<b></HTML>");
        }
    }
    private static String help =
            "      Available commands:\n" +
            "      * means unimplemented command.\n" +
            "cd path                       Change remote directory to 'path'\n" +
            "lcd path                      Change local directory to 'path'\n" +
            "chgrp grp path                Change group of file 'path' to 'grp'\n" +
            "chmod mode path               Change permissions of file 'path' to 'mode'\n" +
            "chown own path                Change owner of file 'path' to 'own'\n" +
            "help                          Display this help text\n" +
            "get remote-path [local-path]  Download file\n" +
            "get-resume remote-path [local-path]  Resume to download file.\n" +
            "get-append remote-path [local-path]  Append remote file to local file\n" +
            "lls [ls-options [path]]       Display local directory listing\n" +
            "ln oldpath newpath            Symlink remote file\n" +
            "*lmkdir path                  Create local directory\n" +
            "lpwd                          Print local working directory\n" +
            "ls [path]                     Display remote directory listing\n" +
            "*lumask umask                 Set local umask to 'umask'\n" +
            "mkdir path                    Create remote directory\n" +
            "put local-path [remote-path]  Upload file\n" +
            "put-resume local-path [remote-path]  Resume to upload file\n" +
            "put-append local-path [remote-path]  Append local file to remote file.\n" +
            "pwd                           Display remote working directory\n" +
            "stat path                     Display info about path\n" +
            "exit                          Quit sftp\n" +
            "quit                          Quit sftp\n" +
            "rename oldpath newpath        Rename remote file\n" +
            "rmdir path                    Remove remote directory\n" +
            "rm path                       Delete remote file\n" +
            "symlink oldpath newpath       Symlink remote file\n" +
            "readlink path                 Check the target of a symbolic link\n" +
            "realpath path                 Canonicalize the path\n" +
            "rekey                         Key re-exchanging\n" +
            "compression level             Packet compression will be enabled\n" +
            "version                       Show SFTP version\n" +
            "?                             Synonym for help";

    public class Get extends Thread {

        String p1;
        String p2;
        SftpProgressMonitor monitor;
        int mode;

        public Get(String p1, String p2, SftpProgressMonitor monitor, int mode) {
            this.p1 = p1;
            this.p2 = p2;
            this.monitor = monitor;
            this.mode = mode;
            this.start();
        }

        @Override
        public void run() {
            try {
                c.get(p1, p2, monitor, mode);
            } catch (SftpException e) {
                String msg = e.getMessage();
               printERR(e.getMessage());
                // java.lang.StringIndexOutOfBoundsException: String index out of range: -1
                /*if (msg.substring(0, msg.indexOf(':')).equals("java.io.FileNotFoundException")) {
                String msg2 = e.getCause().getMessage();
                printERR(msg2);
                } else {
                printERR("Exception not yet supported: " + msg);
                }*/
            }
        }
    }

    public class Put extends Thread {

        String p1;
        String p2;
        SftpProgressMonitor monitor;
        int mode;

        public Put(String p1, String p2, SftpProgressMonitor monitor, int mode) {
            this.p1 = p1;
            this.p2 = p2;
            this.monitor = monitor;
            this.mode = mode;
            this.start();
        }

        @Override
        public void run() {
            try {
                c.put(p1, p2, monitor, mode);
            } catch (SftpException e) {
                String msg = e.getMessage();
                printERR(e.getMessage());
                // java.lang.StringIndexOutOfBoundsException: String index out of range: -1
                /*if (msg.substring(0, msg.indexOf(':')).equals("java.io.FileNotFoundException")) {
                String msg2 = e.getCause().getMessage();
                printERR(msg2);
                } else {
                printERR("Exception not yet supported: " + msg);
                }*/
            }
        }
    }
}
