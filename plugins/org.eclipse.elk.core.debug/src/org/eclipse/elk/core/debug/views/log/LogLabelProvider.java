/*******************************************************************************
 * Copyright (c) 2019 Kiel University and others.
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
package org.eclipse.elk.core.debug.views.log;

import java.time.format.DateTimeFormatter;

import org.eclipse.elk.core.debug.ElkDebugPlugin;
import org.eclipse.elk.core.debug.model.ExecutionInfo;
import org.eclipse.jface.viewers.DelegatingStyledCellLabelProvider.IStyledLabelProvider;
import org.eclipse.jface.viewers.LabelProvider;
import org.eclipse.jface.viewers.StyledString;
import org.eclipse.swt.graphics.Image;

/**
 * Label provider for debug logs.
 */
public class LogLabelProvider extends LabelProvider implements IStyledLabelProvider {

    /** Path for image used for execution infos with logs */
    private static final String LOG_IMAGE_PATH = "/icons/log.png";
    /** Path for image used for execution infos without logs */
    private static final String NO_LOG_IMAGE_PATH = "/icons/no_log.png";

    /** The image used for execution infos with logs. */
    private Image logImage;
    /** The image used for execution infos without logs. */
    private Image noLogImage;

    public LogLabelProvider() {
        // loading icons for different log containers
        logImage =
                ElkDebugPlugin.imageDescriptorFromPlugin(ElkDebugPlugin.PLUGIN_ID, LOG_IMAGE_PATH).createImage();
        noLogImage =
                ElkDebugPlugin.imageDescriptorFromPlugin(ElkDebugPlugin.PLUGIN_ID, NO_LOG_IMAGE_PATH).createImage();
    }

    @Override
    public void dispose() {
        super.dispose();
        
        if (noLogImage != null) {
            noLogImage.dispose();
            noLogImage = null;
        }
        
        if (logImage != null) {
            logImage.dispose();
            logImage = null;
        }
    }

    @Override
    public StyledString getStyledText(Object element) {
        if (element instanceof ExecutionInfo) {
            ExecutionInfo execution = (ExecutionInfo) element;
            StyledString text = new StyledString(execution.getName());
            
            if (execution.getParent() == null) {
                // Add the creation time for top-level elements
                text.append(" (" + DateTimeFormatter.ISO_INSTANT.format(execution.getCreationTime()) + ")",
                        StyledString.COUNTER_STYLER);
            }
            
            return text;
        } else {
            return null;
        }
    }

    @Override
    public Image getImage(Object element) {
        if (element instanceof ExecutionInfo) {
            return ((ExecutionInfo) element).hasLogMessages()
                    ? logImage
                    : noLogImage;
        } else {
            return null;
        }
    }
    
}
