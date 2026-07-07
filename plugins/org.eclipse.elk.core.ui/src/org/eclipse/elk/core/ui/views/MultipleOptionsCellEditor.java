/*******************************************************************************
 * Copyright (c) 2011, 2018 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.ui.views;

import org.eclipse.jface.viewers.ArrayContentProvider;
import org.eclipse.jface.viewers.DialogCellEditor;
import org.eclipse.jface.viewers.LabelProvider;
import org.eclipse.jface.window.Window;
import org.eclipse.swt.events.MouseAdapter;
import org.eclipse.swt.events.MouseEvent;
import org.eclipse.swt.widgets.Composite;
import org.eclipse.swt.widgets.Control;
import org.eclipse.ui.dialogs.ListSelectionDialog;

/**
 * A cell editor that opens a dialog to select a subset from a given number of options.
 */
public class MultipleOptionsCellEditor extends DialogCellEditor {
    /**
     * The items contained in the enumeration. Used to populate the selection dialog.
     */
    private String[] items;
    /**
     * If this is true, the items have been created for an enumeration set. This means that each item is of the
     * following form: {@code EnumConstantName [" " Decorator]}. The optional decorator describes if the enum
     * constant is advanced or experimental and needs to be stripped away, if present, before returning the
     * chosen option values.
     */
    private boolean decoratedItemValues;
    
    
    /**
     * Creates a layouter hint cell editor.
     * 
     * @param parent the parent composite
     * @param items array of items to select from.
     * @param enumerationSet {@code true} if this editor is created for editing the value of an enumeration set.
     */
    public MultipleOptionsCellEditor(final Composite parent, final String[] items, final boolean enumerationSet) {
        super(parent);
        
        this.items = items;
        decoratedItemValues = enumerationSet;
    }
    
    @Override
    protected Control createContents(final Composite cell) {
        Control label = super.createContents(cell);
        label.addMouseListener(new MouseAdapter() {
            public void mouseDoubleClick(final MouseEvent e) {
                Object newValue = openDialogBox(cell);
                if (newValue != null) {
                    markDirty();
                    doSetValue(newValue);
                    fireApplyEditorValue();
                } else {
                    fireCancelEditor();
                }
                // set focus on the layout view in order to be able to respond to key bindings
                LayoutViewPart layoutView = LayoutViewPart.findView();
                if (layoutView != null) {
                    layoutView.setFocus();
                }
            }
        });
        return label;
    }

    @Override
    protected void updateContents(final Object value) {
        if (value instanceof String[]) {
            String[] strArray = (String[]) value;
            
            if (strArray.length == 0) {
                super.updateContents("");
            } else {
                // Build textual representation
                StringBuilder text = new StringBuilder();
                for (String item : strArray) {
                    text.append(", " + item);
                }
                
                // Drop the leading ", " and update the label
                super.updateContents(text.substring(2));
            }
        }
    }

    @Override
    protected Object openDialogBox(final Control cellEditorWindow) {
        ListSelectionDialog dialog = new ListSelectionDialog(
                cellEditorWindow.getShell(),
                items,
                new ArrayContentProvider(),
                new LabelProvider(),
                null);
        dialog.setInitialSelections(matchSelectionToItems((String[]) getValue()));
        
        if (dialog.open() == Window.OK) {
            Object[] result = dialog.getResult();
            String[] stringResult = new String[result.length];
            
            for (int i = 0; i < result.length; i++) {
                stringResult[i] = stripDecorator((String) result[i]);
            }
            
            return stringResult;
        } else {
            return null;
        }
    }
    
    private String stripDecorator(final String text) {
        if (decoratedItemValues) {
            int spaceIndex = text.indexOf(' ');
            
            if (spaceIndex != -1) {
                return text.substring(0, spaceIndex);
            }
        }
        
        return text;
    }
    
    private String[] matchSelectionToItems(final String[] selection) {
        if (decoratedItemValues) {
            String[] selectedItems = new String[selection.length];
            
            for (int i = 0; i < selectedItems.length; i++) {
                // Find the item that starts with the given selection thingy
                for (String candidate : items) {
                    if (candidate.startsWith(selection[i])) {
                        selectedItems[i] = candidate;
                        break;
                    }
                }
            }
            
            return selectedItems;
        } else {
            return selection;
        }
    }

}
