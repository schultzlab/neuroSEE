function neuroSEE_slackNotify( url, slacktext, username )
    SendSlackNotification( url, slacktext, username, ...
       [], [], [], []);    
end