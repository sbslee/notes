# OpenSSH

## Table of Contents

* [Frequently used commands](#Frequently-used-commands)
* [Creating a channel](#Creating-a-channel)

## Frequently used commands <a name="Frequently-used-commands"></a>

To remove all keys belonging to a hostname:

```
ssh-keygen -R host_name
```

## Creating a channel <a name="Creating-a-channel"></a>

First, open your SSH configuration file:

```
vi ~/.ssh/config
```

Next, add the following:

```
Host host_id
    HostName host_name
    User user_name
```
