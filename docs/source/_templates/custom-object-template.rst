{{ name }} {{ objtype[0] | upper }}{{ objtype[1:] }}
{{ '=' * ( (name + objtype) | length + 1 ) }}

.. currentmodule:: {{ module }}

{% if objtype == 'function' %}
.. autofunction:: {{ objname }}
{% endif %}

{% if objtype == 'class' %}
.. autoclass:: {{ objname }}
   :members:
   :show-inheritance:
   :inherited-members:
{% endif %}
