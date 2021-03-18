# OpenSSH

## Table of Contents

* [Frequently used commands](#Frequently-used-commands)
* [Creating a channel with password](#Creating-a-channel-with-password)
* [Creating a channel without password](#Creating-a-channel-without-password)

## Frequently used commands <a name="Frequently-used-commands"></a>

To remove all keys belonging to a host name:

```
ssh-keygen -R host_name
```

## Creating a channel with password <a name="Creating-a-channel-with-password"></a>

First, open your SSH configuration file:

```
vi ~/.ssh/config
```

Next, add the following:

```
Host host_nickname
    HostName host_name
    User user_name
```

You can now access the server:

```
ssh host_nickname
```

## Creating a channel without password <a name="Creating-a-channel-without-password"></a>
