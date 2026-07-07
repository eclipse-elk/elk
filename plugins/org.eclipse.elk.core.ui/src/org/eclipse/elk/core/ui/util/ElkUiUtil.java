/*******************************************************************************
 * Copyright (c) 2009, 2015 Kiel University and others.
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
package org.eclipse.elk.core.ui.util;

import org.eclipse.elk.core.data.LayoutAlgorithmData;
import org.eclipse.elk.core.data.LayoutMetaDataService;
import org.eclipse.elk.core.data.LayoutOptionData;
import org.eclipse.elk.core.options.CoreOptions;
import org.eclipse.emf.common.command.AbstractCommand;
import org.eclipse.emf.edit.domain.EditingDomain;
import org.eclipse.emf.transaction.RecordingCommand;
import org.eclipse.emf.transaction.TransactionalEditingDomain;

/**
 * Utility methods used for the ELK UI.
 *
 * @author msp
 */
public final class ElkUiUtil {
    
    /**
     * Hidden constructor.
     */
    private ElkUiUtil() {
    }
    
    /**
     * Performs the model changes specified in the given runnable in a safe context.
     * If the given editing domain is a {@link TransactionalEditingDomain}, the changes done in
     * the runnable are recorded, otherwise the operation is not undoable.
     * 
     * @param runnable a runnable that performs model changes
     * @param editingDomain the editing domain for the changes, or {@code null}
     * @param label a user friendly label shown for the undo action, or {@code null}
     */
    public static void runModelChange(final Runnable runnable,
            final EditingDomain editingDomain, final String label) {
        if (editingDomain instanceof TransactionalEditingDomain) {
            // execute with a transactional editing domain
            editingDomain.getCommandStack().execute(new RecordingCommand(
                    (TransactionalEditingDomain) editingDomain, label) {
                protected void doExecute() {
                    runnable.run();
                }
            });
        } else if (editingDomain != null) {
            // execute with an arbitrary editing domain
            editingDomain.getCommandStack().execute(new AbstractCommand(label) {
                public void execute() {
                    runnable.run();
                }
                @Override
                public boolean canUndo() {
                    return false;
                }
                public void redo() {
                    execute();
                }
            });
        } else {
            // execute without an editing domain
            runnable.run();
        }
    }
    
    /**
     * Returns the layout option data that matches the given user-friendly name and is known by the
     * given layout provider.
     * 
     * @param providerData a layout provider data
     * @param optionName user-friendly name of a layout option
     * @return the corresponding layout option data
     */
    public static LayoutOptionData getOptionData(final LayoutAlgorithmData providerData,
            final String optionName) {
        for (LayoutOptionData data : LayoutMetaDataService.getInstance().getOptionData()) {
            if (data.getName().equals(optionName) && providerData.knowsOption(data)) {
                return data;
            }
        }
        return null;
    }

    /**
     * Retrieves a suitable layout option data instance that matches the given
     * user friendly display name.
     * 
     * @param providerDataArray array of applicable layout provider data
     * @param displayName display name of the layout option as seen by the user
     * @return the most suitable layout option data
     */
    public static LayoutOptionData getOptionData(final LayoutAlgorithmData[] providerDataArray,
            final String displayName) {
        for (LayoutAlgorithmData providerData : providerDataArray) {
            LayoutOptionData optionData = getOptionData(providerData, displayName);
            if (optionData != null) {
                return optionData;
            }
        }
        // the only option data that is added without explicit support by layouters is layout hint
        return LayoutMetaDataService.getInstance().getOptionData(CoreOptions.ALGORITHM.getId());
    }

}
