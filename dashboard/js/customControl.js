L.Control.Custom = L.Control.Layers.extend({
    onAdd: function () {
        this._controlInputs = [];
        this._initLayout();
        this._base();
        this._addButton();
        this._update();
        this._refresh();
        this._isEnabled = true;
        this._allSelected = false;
        this._noneSelected = false;
        this._checkBoxesCounts = 0;
        return this._container;
    },

    _radioController: function () {
        var selected = this._getSelected();
        if (selected.length === this._checkBoxesCounts) {
            $("#1").prop("checked", true)
            console.log('All checked')
        }

        if (selected.length === 0) {
            $("#0").prop("checked", true)
            console.log('All empty')
        }

        if ([1, this._checkBoxesCounts - 1].includes(selected.length)) {
            $('.control_radio').prop('checked', false);
        }

    },

    _radioMaker: function (id, isChecked) {
        var rb = document.createElement('input');
        rb.type = 'radio';
        rb.className = 'leaflet-control-layers-selector control_radio';
        rb.name = 'control_radio';
        rb.id = id;
        if (isChecked) {
            rb.setAttribute('checked', 'checked')
        }

        return rb
    },

    _base: function () {
        var elements = this._container.getElementsByClassName('leaflet-control-layers-list');
        var baseDiv = L.DomUtil.create('div', 'leaflet-control-layers-base', elements[0]);

        baseDiv.innerHTML = ' <label><div> ' + this._radioMaker('1', true).outerHTML +
            ' <span> Select All  </span>' +
            ' </div></label> ' +
            '<label><div> ' + this._radioMaker('0', false).outerHTML +
            ' <span> Select None  </span>' +
            ' </div></label>' +
            ' <div class="leaflet-control-layers-separator"></div> ';

        // this._controlInputs.push(cb);
        L.DomEvent.on(baseDiv, 'click', this._onRadioClick, this);
    },

    _onRadioClick: function (e) {
        var id = +$("input[name='control_radio']:checked").attr('id');
        id === 1 ? this._checkAll() : this._unCheckAll();

        this._radioController();
        this._refresh();
    },

    _checkAll: function () {
        var inputs = this._getAllInputs();
        for (var i = 0; i < inputs.length; i++) {
            var input = inputs[i];
            if (!input.checked) {
                input.checked = true
            }
        }
    },

    _unCheckAll: function () {
        var inputs = this._getAllInputs();
        for (var i = 0; i < inputs.length; i++) {
            var input = inputs[i];
            if (input.checked) {
                input.checked = false
            }
        }

    },

    _getAllInputs: function () {
        return document.querySelectorAll('input.leaflet-customcontrol-layers-selector' +
            ':not([name="control_radio"])');
    },

    _addButton: function (containerName) {
        if (this._isEnabled) {

            var elements = this._container.getElementsByClassName('leaflet-control-layers-list');
            var div = L.DomUtil.create('div', 'leaflet-control-layers-overlays', elements[0]);
            var cb = document.createElement('input');
            cb.type = 'checkbox';
            cb.className = 'leaflet-customcontrol-layers-selector';
            cb.defaultChecked = 'checked';
            cb.value = containerName;
            cb.name = 'control_checkbox';
            cb.checked = true;

            div.innerHTML = ' <label> ' + cb.outerHTML +
                ' <span> ' + containerName + ' </span>' +
                ' </label> ';

            // this._controlInputs.push(cb);
            L.DomEvent.on(div, 'click', this._onInputClick, this);

            this._checkBoxesCounts = this._checkBoxesCounts + 1

        }
    },

    _contentToggle: function (target) {
        var tempContainer;
        var containerName = target.value;
        if (target.checked) {
            tempContainer = cellContainer_array.filter(d => d.name === containerName)[0];
            // pixiContainer.addChild(tempContainer);
            tempContainer.visible = true
            masterCellRenderer.render(masterCellContainer)
        } else {
            tempContainer = masterCellContainer.getChildByName(containerName);
            tempContainer.visible = false
            // pixiContainer.removeChild(tempContainer);
            masterCellRenderer.render(masterCellContainer)
        }
    },

    _onInputClick: function (e) {
        // var containerName = e.target.value;
        masterCellRenderer.render(masterCellContainer); // Not sure if I need that here. I cant remember if its on purpose or not...it is called a few lines below anyway
        this._contentToggle(e.target);

        this._radioController()

    },

    _refresh: function () {
        //This deserves its own space, it is not a simple function to sit inside the code for this layer control
        var inputs = this._getAllInputs();

        for (var i = 0; i < inputs.length; i++) {
            var input = inputs[i];
            this._contentToggle(input);
        }
    },

    _getSelected: function () {
        var inputs = this._getAllInputs();

        var selected = Array.from(inputs).filter(d => d.checked === true).map(d => d.value);
        return selected
    },

});

L.control.custom = function() {
    return new L.Control.Custom();
}