from libcasm.local_configuration import (
    OccEventSymInfo,
)

from ._clexmonte_system import (
    System,
)


def make_system_event_info(
    system: System,
):
    """Make event info (:class:`libcasm.local_configuration.OccEventSymInfo`) for all
    events in a system

    Parameters
    ----------
    system : libcasm.clexmonte.System

    Returns
    -------
    event_info: dict[str, libcasm.local_configuration.OccEventSymInfo]
        A dictionary of event info
        (:class:`libcasm.local_configuration.OccEventSymInfo`) for each event type in
        the system, with the event type name as the key.
    """
    event_info = dict()
    event_system = system.event_system
    for event_type_name in system.event_type_names:
        (
            phenomenal_clusters,
            equivalent_generating_op_indices,
            translations,
        ) = system.equivalents_info(event_type_name)

        event_info[event_type_name] = OccEventSymInfo.init(
            prim=system.prim,
            system=event_system,
            prototype_event=system.prototype_event(event_type_name),
            phenomenal_clusters=phenomenal_clusters,
            equivalent_generating_op_indices=equivalent_generating_op_indices,
        )

    return event_info
