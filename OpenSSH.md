# OpenSSH

## Table of Contents

* [Frequently used commands](#Frequently-used-commands)
* [Creating a channel with password](#Creating-a-channel-with-password)
* [Creating a channel without password](#Creating-a-channel-without-password)

## Frequently used commands <a name="Frequently-used-commands"></a>

To list all identities represented by the authentication agent:

```
ssh-add -L
```

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
Host host_id
    HostName host_name
    User user_name
```

You can now access the server by:

```
ssh host_id
```

## Creating a channel without password <a name="Creating-a-channel-without-password"></a>

First, run the following:

```
ssh-keygen -t rsa -b 4096 -C "host_id"
```

Save the private key as `host_id_rsa` and the public key as `host_id_rsa.pub`.

Next, add the private key to the authentication agent:

```
ssh-add ~/.ssh/host_id_rsa
```

Add the public key to the server:

```
cat ~/.ssh/host_id_rsa.pub | ssh host_id 'cat >> ~/.ssh/authorized_keys'
```
