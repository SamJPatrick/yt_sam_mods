# yt_thingie

This is a bare-bones template for making a yt extension
package. Minimally, do the following things to set yourself up for
success.

1. Decide on a name for your package. If the name starts with "yt_",
   then you will be able to import it as `yt.extensions.thingie`.
2. Change all the instances of "yt_thingie" to the name of your
   package:
    * the directory with your source code. Currently, it only has a
      file called `__init__.py` in it.
    * In the file `setup.py`, where it says `name=` and `packages=`,
      and also the url to where you have it on github.
    * In the `LICENSE.txt` file.
3. Add your name to `LICENSE.txt` and name/email to `setup.py`.
4. Install it by doing ``pip install -e .`` in the top directory.

## Importing

All source files should go in the inner directory named after your
package. If you have a file in there called `something.py` with a
function called `my_function`. You can import it as:

```
>>> from yt.extensions.thingie.something import my_function
```

If you put the above line into the ``__init__.py`` file. Then you can
also do:
```
>>> from yt.extensions.thingie import my_function
```

The above will only work if you have ``yt`` installed and your package
name starts with "yt_". Otherwise, you can always import your package
as:
```
>>> from PACKAGE.something import my_function
```

## Enjoy!

Just that!
