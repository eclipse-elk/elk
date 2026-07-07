/*******************************************************************************
 * Copyright (c) 2013, 2015 Kiel University and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.

 * This Source Code may also be made available under the following Secondary
 * Licenses when the conditions for such availability set forth in the Eclipse
 * Public License v. 2.0 are satisfied: GPL-3.0 which is available at
 * https://www.gnu.org/licenses/gpl-3.0-standalone.html.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.ui;

import java.util.EnumSet;

import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.elk.core.ui.views.LayoutPropertyDescriptor;
import org.eclipse.jface.resource.ImageRegistry;
import org.eclipse.jface.viewers.LabelProvider;
import org.eclipse.swt.graphics.Image;

/**
 * A label provider for values of layout options. An instance of this label provider must be
 * associated with a specific option in order to correctly translate the passed values.
 *
 * @author msp
 */
public class LayoutOptionLabelProvider extends LabelProvider {
    
    /** the layout option data instance associated with this label provider. */
    private final LayoutOptionData optionData;

    /**
     * Create a label provider for the given layout option.
     * 
     * @param optionData a layout option
     */
    public LayoutOptionLabelProvider(final LayoutOptionData optionData) {
        this.optionData = optionData;
    }
    
    @Override
    public Image getImage(final Object element) {
        ImageRegistry registry = ElkUiPlugin.getInstance().getImageRegistry();
        switch (optionData.getType()) {
        case OBJECT:
        case STRING:
            return registry.get(ElkUiPlugin.IMG_TEXT);
        case BOOLEAN:
            boolean istrue = true;
            if (element instanceof Boolean) {
                istrue = (Boolean) element;
            } else if (element instanceof Integer) {
                istrue = (Integer) element == 1;
            }
            if (istrue) {
                return registry.get(ElkUiPlugin.IMG_TRUE);
            } else {
                return registry.get(ElkUiPlugin.IMG_FALSE);
            }
        case ENUM:
        case ENUMSET:
            return registry.get(ElkUiPlugin.IMG_CHOICE);
        case INT:
            return registry.get(ElkUiPlugin.IMG_INT);
        case DOUBLE:
            return registry.get(ElkUiPlugin.IMG_DOUBLE);
        default:
            return null;
        }
    }

    @Override
    @SuppressWarnings("rawtypes")
    public String getText(final Object element) {
        switch (optionData.getType()) {
        case STRING:
            if (CoreOptions.ALGORITHM.equals(optionData) 
                    || optionData.getId().equals(LayoutPropertyDescriptor.DISCO_LAYOUT_ALG_ID)) {
                LayoutMetaDataService layoutDataService = LayoutMetaDataService.getInstance();
                LayoutAlgorithmData algorithmData = layoutDataService.getAlgorithmData((String) element);
                if (algorithmData != null) {
                    String bundleName = algorithmData.getBundleName();
                    if (bundleName == null) {
                        return algorithmData.getName();
                    } else {
                        return algorithmData.getName() + " (" + bundleName + ")";
                    }
                }
                return Messages.getString("elk.ui.8");
            }
            break;
        case BOOLEAN:
            if (element instanceof Boolean) {
                return ((Boolean) element).toString();
            }
            // fall through so the same method as for enums is applied
        case ENUM:
            if (element instanceof Integer) {
                return optionData.getChoices()[(Integer) element];
            }
            break;
        case ENUMSET:
            if (element instanceof String) {
                return (String) element;
            } else if (element instanceof String[]) {
                String[] arr = (String[]) element;
                if (arr.length == 0) {
                    return "";
                } else {
                    StringBuilder builder = new StringBuilder();
                    
                    for (String s : arr) {
                        builder.append(", ").append(s);
                    }
                    
                    return builder.substring(2);
                }
            } else if (element instanceof EnumSet) {
                EnumSet set = (EnumSet) element;
                if (set.isEmpty()) {
                    return "";
                }
                
                StringBuilder builder = new StringBuilder();
                for (Object o : set) {
                    builder.append(", " + ((Enum) o).name());
                }
                return builder.substring(2);
            }
        }
        return element.toString();
    }
    
}
