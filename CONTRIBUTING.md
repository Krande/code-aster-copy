# To contribute to the project, you need to do it through a merge request

Before making your first contributions, please sign and sent the contribution
agreement available from [this page][1].

First you need to fork the repository into your own account. You can do that
simply by clicking the **fork** button on the GitLab interface.

Then, clone the repository on your laptop.

```shell
git clone https://gitlab.com/your-username/your-forkname.git
```

Once this is done, you can setup the *codeaster/src* repository as the *upstream*
of your clone to simplify the update of your fork repository.

```shell
git remote add upstream https://gitlab.com/codeaster/src.git
```

Now, you have your repository configured, and you want to create a new merge request.
The first step is to create a branch from the HEAD of the **main** branch of
your fork repository.

```shell
git checkout -b your_branch_name
```

Apply your modifications in your branch. Then, you need to push this branch on
your online repository.

```shell
git push origin your_branch_name
```

or, if the branch already exists online, and you want to update it:

```shell
git push -f origin your_branch_name
```

Once your branch is online, on the GitLab interface, go to the branches webpage,
select the branch you want to push as a merge request, and push the button !

***Be careful to check the 'close after merge' check box, and to push to the
codeaster/src repository***.

[1]: https://gitlab.com/codeaster-opensource-documentation/opensource-installation-development/-/blob/main/devel/merge_request.md
