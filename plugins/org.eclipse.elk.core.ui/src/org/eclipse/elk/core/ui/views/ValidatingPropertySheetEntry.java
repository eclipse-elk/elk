/*******************************************************************************
 * Copyright (c) 2016 TypeFox and others.
 * 
 * This program and the accompanying materials are made available under the
 * terms of the Eclipse Public License 2.0 which is available at
 * http://www.eclipse.org/legal/epl-2.0.
 *
 * SPDX-License-Identifier: EPL-2.0 OR GPL-3.0-or-later
 *******************************************************************************/
package org.eclipse.elk.core.ui.views;

import java.util.Collections;
import java.util.List;

import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.core.ui.ElkUiPlugin;
import org.eclipse.elk.core.validation.GraphIssue;
import org.eclipse.jface.viewers.DecorationOverlayIcon;
import org.eclipse.jface.viewers.IDecoration;
import org.eclipse.swt.graphics.Image;
import org.eclipse.ui.views.properties.IPropertyDescriptor;
import org.eclipse.ui.views.properties.IPropertySource;
import org.eclipse.ui.views.properties.PropertySheetEntry;

/**
 * A property sheet entry that shows error markers and error messages when validation issues
 * are found for a layout option value.
 */
public class ValidatingPropertySheetEntry extends PropertySheetEntry {
    
    @Override
    protected PropertySheetEntry createChildEntry() {
        return new ValidatingPropertySheetEntry();
    }
    
    /**
     * @return the property source
     */
    protected LayoutPropertySource getSource() {
        ValidatingPropertySheetEntry e = this;
        do {
            Object[] values = e.getValues();
            if (values.length > 0) {
                IPropertySource source = e.getPropertySource(values[0]);
                if (source instanceof LayoutPropertySource) {
                    return (LayoutPropertySource) source;
                }
            }
            e = (ValidatingPropertySheetEntry) e.getParent();
        } while (e != null);
        return null;
    }
    
    /**
     * Validate the displayed value by delegating to the {@link LayoutPropertySource}, if any.
     */
    protected List<GraphIssue> validate() {
        LayoutPropertySource source = getSource();
        IPropertyDescriptor descriptor = getDescriptor();
        Object[] values = getValues();
        if (source != null && descriptor instanceof LayoutPropertyDescriptor && values.length > 0) {
            LayoutOptionData optionData = ((LayoutPropertyDescriptor) descriptor).getOptionData();
            Object optionValue = source.translateFromUI(values[0], optionData);
            List<GraphIssue> result = source.validatePropertyValue(optionData, optionValue);
            if (result != null) {
                return result;
            }
        }
        return Collections.emptyList();
    }
    
    @Override
    public Image getImage() {
        Image image = super.getImage();
        List<GraphIssue> issues = validate();
        if (issues.isEmpty()) {
            return image;
        } else if (image != null) {
            String overlayKey = "error_overlay/" + image.toString();
            Image result = ElkUiPlugin.getInstance().getImageRegistry().get(overlayKey);
            if (result == null) {
                DecorationOverlayIcon overlayIcon = new DecorationOverlayIcon(image,
                        ElkUiPlugin.getImageDescriptor("icons/error_overlay.png"), IDecoration.BOTTOM_LEFT);
                result = overlayIcon.createImage();
                ElkUiPlugin.getInstance().getImageRegistry().put(overlayKey, result);
            }
            return result;
        }
        return null;
    }
    
    @Override
    public String getDescription() {
        String description =  super.getDescription();
        List<GraphIssue> issues = validate();
        if (issues.isEmpty()) {
            return description;
        } else {
            StringBuilder builder = new StringBuilder();
            GraphIssue firstIssue = issues.get(0);
            builder.append(firstIssue.getSeverity().getUserString())
                .append(": ")
                .append(firstIssue.getMessage());
            if (issues.size() > 1) {
                builder.append(" (").append(issues.size() - 1).append(" more issues)");
            }
            return builder.toString();
        }
    }

}
