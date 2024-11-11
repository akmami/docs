# SSHFS README

## Overview

SSHFS (SSH Filesystem) allows you to mount a remote filesystem on your local machine over SSH. This makes it possible to work with files on a remote server as though they are on your local system, all with secure SSH encryption.

SSHFS is particularly useful for developers and system administrators who need to manage files on remote servers seamlessly.

## Prerequisites

- SSH access to the remote server
- SSHFS installed on your local machine

## Installation

### On Ubuntu/Debian

```bash
sudo apt update
sudo apt install sshfs
```

### On macOS (using Homebrew)

```bash
brew install sshfs
```

## Usage

### Basic Command Syntax

```bash
sshfs [user]@[host]:[remote_directory] [local_mount_point]
```

- `[user]`: Username on the remote server
- `[host]`: IP address or hostname of the remote server
- `[remote_directory]`: Directory on the remote server you want to mount
- `[local_mount_point]`: Directory on your local machine where the remote directory will be accessible

### Example Commands

1. **Mount a Remote Directory**

   ```bash
   sshfs user@remote_host:/path/to/remote/directory /path/to/local/mountpoint
   ```

2. **Unmount the Mounted Directory**

   ```bash
   fusermount -u /path/to/local/mountpoint   # Linux
   umount /path/to/local/mountpoint          # macOS
   ```

### Mount with Specific Options

You can add options using the `-o` flag. Some useful options include:

- **reconnect**: Automatically reconnect if the connection drops.
- **IdentityFile**: Specify an SSH key file.
- **allow_other**: Allow other users to access the mounted directory (requires setting user_allow_other in `/etc/fuse.conf`).

**Example with Options:**

```bash
sshfs -o IdentityFile=~/.ssh/id_rsa,reconnect user@remote_host:/path/to/remote /path/to/local_mount
```

## Tips

- **Permissions**: Make sure the user has the necessary permissions on the remote directory.
- **SSH Keys**: If possible, use SSH key authentication for a smoother experience.
