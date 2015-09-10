% Email settings
setpref('Internet','E_mail','o.komarov11@imperial.ac.uk');
setpref('Internet','SMTP_Server','smtp.cc.ic.ac.uk')
setpref('Internet','SMTP_Username','ok1011');
if ~ispref('Internet','SMTP_Password')
    pwd = inputdlg('Email pwd');
    setpref('Internet','SMTP_Password',pwd{:});
end
props = java.lang.System.getProperties;
props.setProperty('mail.smtp.auth','true');
props.setProperty('mail.smtp.socketFactory.class', ...
                  'javax.net.ssl.SSLSocketFactory');
props.setProperty('mail.smtp.socketFactory.port','465');
clear props pwd